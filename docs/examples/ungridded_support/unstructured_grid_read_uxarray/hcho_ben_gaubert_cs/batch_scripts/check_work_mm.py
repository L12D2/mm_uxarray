
"""Audit MELODIES-MONET changes

Usage (glade):
    python audit_paired_output.py '<output_dir>/202406*_tropomi_l2_co_cam-chem-se*.nc4'
    python audit_paired_output.py '<output_dir>/202406*.nc4' --days 30

Exit code 1 if any FAIL, so it can gate a plotting job:
    python audit_paired_output.py "$OUT/202406*.nc4" && qsub plot_sat.sh
"""
import argparse
import glob
import os
import sys

import numpy as np
import xarray as xr

# plausible physical ranges for domain-median columns, molec/cm2
MEDIAN_FENCES = {
    "vertical_column":                            (1e14, 5e17),  # TEMPO HCHO
    "formaldehyde_tropospheric_vertical_column":  (1e14, 5e17),  # TROPOMI HCHO
    "vertical_column_troposphere":                (1e13, 5e17),  # TEMPO NO2
    "nitrogendioxide_tropospheric_column":        (1e13, 5e17),  # TROPOMI NO2
    "carbonmonoxide_total_column":                (5e17, 1e19),  # TROPOMI CO
    "CH2O": (1e14, 5e17), "NO2": (1e13, 5e17), "CO": (5e17, 1e19),
    # surface, ppb
    "O3": (0.1, 300), "O3_new": (0.1, 300),
    "NO2_new": (0.001, 500), "PM2.5": (0.01, 500),
}

SFC_FENCES = {
    "O3": (0.1, 300), "NO2": (0.01, 500), "CO": (0.005, 5000),
    "PM2.5": (0.01, 500), "TEMP": (200, 330),
}

# obs var to model var expected in the same pair file
OBS2MOD = {
    "vertical_column": "CH2O",
    "formaldehyde_tropospheric_vertical_column": "CH2O",
    "vertical_column_troposphere": "NO2",
    "nitrogendioxide_tropospheric_column": "NO2",
    "carbonmonoxide_total_column": "CO",
}

CONUS = (-130.0, -60.0, 20.0, 55.0)

def _fail(msgs, tag, text):
    msgs.append(("FAIL", tag, text))

def _warn(msgs, tag, text):
    msgs.append(("WARN", tag, text))

def audit_file(path):
    msgs = []
    try:
        ds = xr.open_dataset(path)
    except Exception as e:
        _fail(msgs, "open", f"unreadable (corrupt/killed write?): {e!r}")
        return msgs
    with ds:
        is_obsgrid = "_obsgrid" in path
        is_series = "_series" in path

        # --- grid shape / extent regressions -----------------------------
        if is_obsgrid:
            lon = ds["longitude"].values if "longitude" in ds else None
            if lon is None:
                _fail(msgs, "grid", "obsgrid file without longitude")
            else:
                lo, hi = float(np.nanmin(lon)), float(np.nanmax(lon))
                if hi - lo > 180:
                    _fail(msgs, "grid",
                          f"GLOBAL grid (lon {lo:.1f}..{hi:.1f}) — "
                          "extent fix missing when this was paired; re-pair")
                elif not (CONUS[0] - 5 <= lo and hi <= CONUS[1] + 5):
                    _warn(msgs, "grid",
                          f"lon {lo:.1f}..{hi:.1f} outside CONUS+5deg")

        # --- variable presence / units / medians --------------------------
        seen_pairable = False
        for ov, mv in OBS2MOD.items():
            if ov not in ds:
                continue
            seen_pairable = True
            if mv not in ds:
                _fail(msgs, "vars", f"obs {ov} present but model {mv} missing")
                continue
            o = ds[ov].values
            m = ds[mv].values
            o_med = float(np.nanmedian(o))
            m_med = float(np.nanmedian(m))
            for name, med in ((ov, o_med), (mv, m_med)):
                fence = MEDIAN_FENCES.get(name)
                if fence and np.isfinite(med) and not (fence[0] <= abs(med) <= fence[1]):
                    _fail(msgs, "units",
                          f"{name} median {med:.3e} outside {fence} — "
                          "units/vertical-coordinate slip?")
            if np.isfinite(o_med) and np.isfinite(m_med) and o_med != 0:
                r = m_med / o_med
                if not (0.05 <= abs(r) <= 20):
                    _fail(msgs, "ratio",
                          f"model/obs median ratio {r:.2f} ({mv}/{ov}) — "
                          "operator or units problem")
                elif not (0.2 <= abs(r) <= 5):
                    _warn(msgs, "ratio",
                          f"model/obs median ratio {r:.2f} ({mv}/{ov})")
            units = ds[ov].attrs.get("units", "")
            if is_obsgrid or is_series:
                if units and "cm" not in units.replace(" ", ""):
                    _warn(msgs, "units",
                          f"{ov} units attr {units!r} (expected molec/cm2)")
            # all-NaN / coverage
            finite = np.isfinite(o).mean() if o.size else 0.0
            if finite == 0:
                _fail(msgs, "data", f"{ov} is all-NaN")
            elif finite < 0.005 and not is_series:
                _warn(msgs, "data", f"{ov} only {finite:.2%} finite")

        # surface pairs (no satellite obs var matched)
        if not seen_pairable:
            dvars = set(ds.data_vars)
            sfc = [v for v in dvars if v in SFC_FENCES]
            for v in sfc:
                o_med = float(np.nanmedian(ds[v].values))
                fence = SFC_FENCES[v]
                if np.isfinite(o_med) and not (fence[0] <= abs(o_med) <= fence[1]):
                    _warn(msgs, "range",
                          f"{v} median {o_med:.3g} outside {fence}")
                mv = v + "_new" if v + "_new" in dvars else None
                if mv:
                    m_med = float(np.nanmedian(ds[mv].values))
                    if np.isfinite(o_med) and np.isfinite(m_med) and o_med != 0:
                        r = m_med / o_med
                        if not (0.05 <= abs(r) <= 20):
                            _fail(msgs, "units",
                                  f"model/obs median ratio {r:.1f} "
                                  f"({mv} {m_med:.3g} / {v} {o_med:.3g}) — "
                                  "ppm/ppb conversion missing?")
            if not sfc:
                _warn(msgs, "vars",
                      f"no recognized variables; has {sorted(dvars)[:8]}")

        # --- series files need obs error for inverse-variance weighting ---
        if is_series:
            for v in ds.data_vars:
                w = ds[v].attrs.get("series_weighting")
                if w is None:
                    _warn(msgs, "series",
                          f"{v}: no series_weighting attr (paired before the "
                          "attr existed) — weighting unverifiable; re-pair to "
                          "stamp it")
                    break
                if "fallback" in w:
                    _warn(msgs, "series",
                          f"{v}: {w} — add <obs_var>_precision/_uncertainty "
                          "to the obs variables in the control")

        # --- time sanity ---------------------------------------------------
        if "time" in ds.coords and ds["time"].size:
            t = np.asarray(ds["time"].values)
            if np.unique(t).size != t.size:
                _warn(msgs, "time", "duplicate time stamps (double-paired day?)")
    return msgs

def main():
    p = argparse.ArgumentParser()
    p.add_argument("patterns", nargs="+", help="globs of paired .nc4 files")
    p.add_argument("--days", type=int, default=None,
                   help="expected number of distinct YYYYMMDD prefixes "
                        "(e.g. 30 for June); checks for missing days")
    args = p.parse_args()

    files = sorted(set(f for pat in args.patterns for f in glob.glob(pat)))
    if not files:
        print("no files matched", file=sys.stderr)
        sys.exit(2)

    n_fail = 0
    for f in files:
        msgs = audit_file(f)
        fails = [m for m in msgs if m[0] == "FAIL"]
        n_fail += len(fails)
        status = "FAIL" if fails else ("warn" if msgs else "ok")
        print(f"[{status:4s}] {os.path.basename(f)}")
        for lvl, tag, text in msgs:
            print(f"       {lvl} {tag}: {text}")

    # missing-day check per product (filenames start YYYYMMDD_)
    if args.days:
        from collections import defaultdict
        by_prod = defaultdict(set)
        for f in files:
            base = os.path.basename(f)
            day, _, rest = base.partition("_")
            if len(day) == 8 and day.isdigit():
                by_prod[rest].add(day)
        for prod, days in sorted(by_prod.items()):
            if len(days) < args.days:
                have = sorted(days)
                print(f"[warn] {prod}: {len(days)}/{args.days} days "
                      f"({have[0]}..{have[-1]})")
                n_fail += 0  # coverage gaps warn, don't fail

    print(f"\n{len(files)} files, {n_fail} FAILs")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()