"""
Need a way to compare different obs simultaneously rather than one at a time. Because of how involved this may become, i wrote a septerate 
obs2obs util folder to keep any issues with it seperate (as much as possible) from the rest of the MM 


obs2obs: compare multiple obs AND multiple models simultaneously.

MELODIES-MONET pairing is one-obs to N-models; this util works AFTER pairing,
on the saved pair objects (an.paired), treating each (pair label, variable) as
an independent *data series*. Because the satellite obs-grid pairs share one
lat/lon grid (and model-space pairs share the model mesh)

multi-obs comparison should be pure alignment -- no new pairing needed 

Series spec (each entry in a group's ``series`` list):
    - {label: <pair label>, var: <variable>, name: <legend>, color: <opt>}
      raw values ("value" mode). Only mix series with the SAME units.
    - {label: ..., obs_var: X, mod_var: Y, name: ..., mode: bias}
      model - obs from within one pair.
    - {label: ..., obs_var: X, mod_var: Y, name: ..., mode: relative_bias}
      100*(model-obs)/obs. Dimensionless units will allow to put surface ppb,
      satellite columns, and aircraft data on one axis.

Handles pair objects that are xarray Datasets (satellite: (time,y,x),
(time,n_face), or (time,) series) and pandas DataFrames (surface/aircraft
point obs). Time matching: 'none', 'daily', or {window_local: [h0, h1]}
(solar local hour ~ UTC + lon/15; the honest way to co-locate a
geostationary full-day product with a ~13:30 LT polar orbiter).

Driver file: 

    an.read_control(); an.read_analysis()
    from melodies_monet.util import obs2obs_util
    obs2obs_util.run(an.paired, an.control_dict["obs2obs"],
                     default_outdir=an.output_dir)
                     
"""

import os
import gc

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# use native plotting capabilities within MM 
from melodies_monet.plots import surfplots as splots
from melodies_monet.plots import xarray_plots as xrplots

from melodies_monet.plots import savefig as _mm_savefig

def _save(outname_png):
    """Native MM save (adds the logo), same call pattern as xarray_plots."""
    try:
        _mm_savefig(outname_png, logo_height=100, bbox_inches="tight", dpi=200)
    except Exception as e:  # never lose a figure over logo trouble
        print(f"obs2obs: mm savefig failed ({e!r}); plain matplotlib save.",
              flush=True)
        plt.savefig(outname_png, dpi=200, bbox_inches="tight")

def _fnum(x):
    """PyYAML parses '5.0e16' (no '+') as a STRING, which then compares
    silently-wrong; coerce every user-supplied numeric defensively."""
    return None if x is None else float(x)

def _clip2(clip):
    """clip: [lo, hi] from YAML (float, float) or None.
       Add an additioanl way for user to clip values they may think are extraneous 
    """
    return None if clip is None else (float(clip[0]), float(clip[1]))
    
def _pair_obj(paired, label):
    if label not in paired:
        raise KeyError(
            f"obs2obs: pair label {label!r} not in paired data. "
            f"Available: {sorted(paired)}. Check read.paired.filenames.")
    return paired[label].obj

def _as_dataframe(obj):
    """Surface/aircraft pairs may carry a (Multi)Index; normalize to columns."""
    return obj.reset_index() if not isinstance(obj.index, pd.RangeIndex) else obj

def _local_hour_mask_xr(da, dset, h0, h1):
    """(time, ...) mask where solar local hour (UTC + lon/15) is in [h0, h1]."""
    if "longitude" not in dset.coords and "longitude" not in dset.variables:
        print("obs2obs: window_local requested but no longitude on this pair; "
              "skipping the window for it.", flush=True)
        return None
    hours = da["time"].dt.hour + da["time"].dt.minute / 60.0
    local = (hours + dset["longitude"] / 15.0) % 24.0
    mask = (local >= h0) & (local <= h1)
    
    try:
        frac = float(mask.mean().values)
        mean_kept = float(local.where(mask).mean().values)
        print(f"obs2obs: window_local [{h0},{h1}] LST -> kept {frac * 100:.0f}% "
              f"of (time x cell) samples, mean kept local hour {mean_kept:.1f} "
              "(TROPOMI should read ~13.5)", flush=True)
    except Exception:  # diagnostics must never sink a plot
        pass
    return mask

def _subsample(flat, max_points):
    """Deterministic strided subsample so month-scale swaths (~1e8 points)
    don't blow up the boxplot frame; quartiles are unaffected."""
    if max_points and flat.size > max_points:
        step = int(np.ceil(flat.size / max_points))
        return flat[::step]
    return flat

def _extract(spec, paired, time_match, max_points=2_000_000):
    """One series: dict(name, color, ts: pandas Series over time,
    flat: 1-D valid values, grid: time-mean (y,x) DataArray or None)."""
    label = spec["label"]
    mode = spec.get("mode", "value")
    name = spec.get("name", f"{label}:{spec.get('var', mode)}")
    obj = _pair_obj(paired, label)

    if isinstance(obj, xr.Dataset):
        if mode == "value":
            da = obj[spec["var"]]
        else:
            da_o, da_m = obj[spec["obs_var"]], obj[spec["mod_var"]]
            da = (da_m - da_o) if mode == "bias" \
                else 100.0 * (da_m - da_o) / da_o.where(da_o > 0)

        # temporary filter string
        clip = _clip2(spec.get("clip"))   # [lo, hi] physical range -> NaN outside
        if clip is not None:
            da = da.where((da >= clip[0]) & (da <= clip[1]))
            
        bounds = spec.get("bounds")
        if bounds is not None and "longitude" in obj.variables:
            w, e, s, n = bounds
            da = da.where((obj["longitude"] >= w) & (obj["longitude"] <= e)
                          & (obj["latitude"] >= s) & (obj["latitude"] <= n))

        if isinstance(time_match, dict) and "window_local" in time_match:
            h0, h1 = time_match["window_local"]
            m = _local_hour_mask_xr(da, obj, h0, h1)
            if m is not None:
                da = da.where(m)

        spatial = [d for d in da.dims if d != "time"]
        ts = da.mean(dim=spatial, skipna=True).to_series() if "time" in da.dims \
            else pd.Series(dtype=float)
        if time_match == "daily" and not ts.empty:
            ts = ts.resample("1D").mean()

        vals = da.values.ravel()
        flat = _subsample(vals[np.isfinite(vals)], max_points)

        grid = None
        if {"y", "x"} <= set(da.dims):
            grid = da.mean(dim="time", skipna=True) if "time" in da.dims else da

        # da is the masked / windowed array itself for the plots that need 
        # point level alignmet between timeseries (e.g. the scatter plots) rather than just a pre-reduced ts/flat/grid view
        return dict(name=name, color=spec.get("color"), ts=ts, flat=flat,
                    grid=grid, da = da, src=obj)

    # pandas path (surface / aircraft point pairs)
    df = _as_dataframe(obj)
    if mode == "value":
        v = df[spec["var"]]
    else:
        o, m = df[spec["obs_var"]], df[spec["mod_var"]]
        v = (m - o) if mode == "bias" else 100.0 * (m - o) / o.where(o > 0)
    
    clip = _clip2(spec.get("clip"))   # [lo, hi] physical range -> NaN outside
    if clip is not None:
        v = v.where((v >= clip[0]) & (v <= clip[1]))
        
    work = pd.DataFrame({"time": pd.to_datetime(df["time"]), "v": v})
    if "longitude" in df.columns:
        work["longitude"] = df["longitude"].values
        work["latitude"] = df["latitude"].values

    bounds = spec.get("bounds")
    if bounds is not None and "longitude" in work.columns:
        w, e, s, n = bounds
        work = work[(work.longitude >= w) & (work.longitude <= e)
                    & (work.latitude >= s) & (work.latitude <= n)]

    if isinstance(time_match, dict) and "window_local" in time_match:
        if "longitude" in work.columns:
            h0, h1 = time_match["window_local"]
            lh = (work.time.dt.hour + work.time.dt.minute / 60.0
                  + work.longitude / 15.0) % 24.0
            work = work[(lh >= h0) & (lh <= h1)]
        else:
            print(f"obs2obs: window_local skipped for {label} (no longitude).",
                  flush=True)

    ts = work.set_index("time")["v"].resample(
        "1D" if time_match == "daily" else "1h").mean().dropna()
    flat = work["v"].values
    flat = _subsample(flat[np.isfinite(flat)], max_points)
    return dict(name=name, color=spec.get("color"), ts=ts, flat=flat,
                grid=None, da=None, src=obj)


# MM plotting 

def _ts_to_dataset(ts):
    """Aligned pandas Series: 1-var Dataset xrplots.make_timeseries accepts."""
    idx = pd.DatetimeIndex(ts.index, name="time")
    da = xr.DataArray(ts.values, dims=("time",), coords={"time": idx}, name="v")
    return da.to_dataset()


def multi_timeseries(series, outname, ylabel=None, title=None, avg_window=None,
                     fig_dict=None, text_dict=None):
    """Overlay every series on one axis via xrplots.make_timeseries."""

    tk = {"fontsize": 16, **(text_dict or {})}
    fig, ax = plt.subplots(figsize=(fig_dict or {}).get("figsize", (12, 6)))
    lo, hi, drawn = np.inf, -np.inf, 0
    
    for s in series:
        ts = s["ts"]
        if avg_window and not ts.empty:
            ts = ts.resample(avg_window).mean()
        ts = ts.dropna()
        if ts.notna().sum() == 0:
            print(f"obs2obs: no timeseries for {s['name']}; skipping.", flush=True)
            continue
        ax.plot(ts.index, ts.values, marker="*", linewidth=2.0,
                color=s["color"], label=s["name"])
        lo = min(lo, float(np.nanmin(ts.values)))
        hi = max(hi, float(np.nanmax(ts.values)))
        drawn += 1
    if not drawn:
        plt.close(fig)
        return
    pad = 0.05 * (hi - lo or 1.0)
    ax.set_ylim(lo - pad, hi + pad)           # span all series, not just obs
    if ylabel:
        ax.set_ylabel(ylabel, fontweight="bold", fontsize=tk["fontsize"])
    ax.set_xlabel("time", fontweight="bold", fontsize=tk["fontsize"])
    
    ax.legend(fontsize=tk["fontsize"] * 0.7)
    if title:
        ax.set_title(title, fontweight="bold", fontsize=tk["fontsize"])
    plt.tight_layout()
    _save(outname + ".png")
    plt.close(fig)
    print(f"obs2obs: wrote {outname}.png ({drawn} series)", flush=True)


def multi_boxplot(series, outname, ylabel=None, title=None, vmin=None,
                  vmax=None, fig_dict=None, text_dict=None, showfliers=False,
                        hline=None): # add outliers as an option in the yaml
    """One box per series via surfplots.make_boxplot.

    The combined frame is built with pd.concat (NOT column assignment as in
    calculate_boxplot) because our series have different lengths.
    """
    cols, label_bx = {}, []
    for s in series:
        if s["flat"].size == 0:
            print(f"obs2obs: no valid data for {s['name']}; skipping box.",
                  flush=True)
            continue
        cols[s["name"]] = pd.Series(s["flat"])
        label_bx.append(dict(color=s["color"] or "gray", linestyle="-",
                             marker="x", linewidth=1.2, markersize=6.0,
                             column=s["name"], label=s["name"]))
    if not cols:
        return
    comb_bx = pd.concat(cols, axis=1)   # union index; shorter series NaN-pad
    splots.make_boxplot(comb_bx, label_bx, ylabel=ylabel, vmin=vmin, vmax=vmax,
                        outname=outname, domain_type="all",
                        domain_name=title or "", fig_dict=fig_dict,
                        text_dict=text_dict, showfliers=showfliers, hline = hline, debug=False)
    print(f"obs2obs: wrote {outname} (boxplot)", flush=True)


def diff_map(series, outname, ylabel=None, title=None, vdiff=None,
             domain_name=None, fig_dict=None, text_dict=None):
    """series[0] - series[1] mean map via xrplots.make_spatial_bias_gridded."""
    a, b = series[0], series[1]
    if a["grid"] is None or b["grid"] is None:
        print("obs2obs: diff_map needs two gridded (y,x) series "
              f"({a['name']}, {b['name']}); skipping.", flush=True)
        return
    if a["grid"].shape != b["grid"].shape:
        print(f"obs2obs: diff_map grids differ {a['grid'].shape} vs "
              f"{b['grid'].shape}; regrid to a common obs grid first.",
              flush=True)
        return
    ds = xr.Dataset(
        {"s_minuend": a["grid"], "s_subtrahend": b["grid"]},
        coords={"longitude": a["src"]["longitude"],
                "latitude": a["src"]["latitude"]},
    )

    if "time" in a["src"].variables or "time" in a["src"].coords:
        _t = np.asarray(a["src"]["time"].values)
        if _t.size:
            ds = ds.assign_coords(time=("time", [_t.min(), _t.max()]))
            
    # plots s_minuend - s_subtrahend; the auto-region branch zooms to the
    # finite footprint of the difference.
    xrplots.make_spatial_bias_gridded(
        ds, varname_o="s_subtrahend", label_o=b["name"],
        varname_m="s_minuend", label_m=a["name"],
        ylabel=ylabel or f"{a['name']} - {b['name']}", vdiff=vdiff,
        outname=outname, domain_type="auto-region:box",
        domain_name=domain_name or (title or "obs2obs"),
        fig_dict=fig_dict or {"states": True, "counties": False},
        text_dict=text_dict, debug=False)
    print(f"obs2obs: wrote {outname} (diff_map)", flush=True)

def scatter(series, outname, ylabel=None, xlabel=None, title=None,
            fig_dict=None, text_dict=None, density=False, xylim=None,
            time_avg="1D", max_plot_points=200_000):
    """Point-matched scatter of series[1] (y) against series[0] (x).

    Alignment before plotting

    - two GRIDDED series (obsgrid pairs): daily-mean per cell (``time_avg``),
      then inner-join on the shared (time, y, x); every point is the same
      cell on the same day seen by both series. Requires the same grid --
      e.g. obs vs model within one pair (always same grid), or two pairs on
      the same obsgrid resolution/extent.
    - otherwise: daily domain-mean timeseries join (one point per day).

    Draws the 1:1 line and an N/r/slope/MB/RMSE stats box (computed on ALL
    matched points), then renders via surfplots.make_scatter_density_plot
    (best-fit line included; ``density: true`` -> seaborn KDE fill).
    
    """
    a, b = series[0], series[1]

    def _daily(s):
        da = s.get("da")
        if isinstance(da, xr.DataArray) and "time" in da.dims and time_avg:
            return da.resample(time=time_avg).mean(skipna=True)
        return da

    xa, ya = _daily(a), _daily(b)
    if isinstance(xa, xr.DataArray) and isinstance(ya, xr.DataArray):
        if (xa.sizes.get("y"), xa.sizes.get("x")) != (ya.sizes.get("y"), ya.sizes.get("x")):
            print(f"obs2obs: scatter {a['name']!r} vs {b['name']!r}: grids "
                  f"differ ({dict(xa.sizes)} vs {dict(ya.sizes)}). Cell-matched "
                  "scatter needs the same obsgrid -- pair both products at the "
                  "same res/extent, or compare within one pair (obs vs model).",
                  flush=True)
            return
        xa, ya = xr.align(xa, ya, join="inner")
        x, y = xa.values.ravel(), ya.values.ravel()
    else:
        ja = a["ts"].resample("1D").mean() if not a["ts"].empty else a["ts"]
        jb = b["ts"].resample("1D").mean() if not b["ts"].empty else b["ts"]
        j = pd.concat([ja, jb], axis=1, join="inner").dropna()
        if j.empty:
            print(f"obs2obs: scatter {a['name']!r} vs {b['name']!r}: no "
                  "overlapping days.", flush=True)
            return
        x, y = j.iloc[:, 0].to_numpy(float), j.iloc[:, 1].to_numpy(float)

    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 3:
        print(f"obs2obs: scatter {a['name']!r} vs {b['name']!r}: only "
              f"{x.size} matched points; skipping.", flush=True)
        return

    # stats on the FULL matched set (plot may be subsampled below)
    n = int(x.size)
    mb = float(np.mean(y - x))
    rmse = float(np.sqrt(np.mean((y - x) ** 2)))
    r = float(np.corrcoef(x, y)[0, 1])
    slope, intercept = (float(v) for v in np.polyfit(x, y, 1))
    xp, yp = _subsample(x, max_plot_points), _subsample(y, max_plot_points)

    if xylim is not None:
        lo, hi = float(xylim[0]), float(xylim[1])
    else:
        lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
        pad = 0.05 * (hi - lo or 1.0)
        lo, hi = lo - pad, hi + pad

    tk = {"fontsize": 14, **(text_dict or {})}
    fig, ax = plt.subplots(figsize=(fig_dict or {}).get("figsize", (7.5, 7)))
    ax.plot([lo, hi], [lo, hi], color="grey", lw=1.2, label="1:1")
    ax.text(0.97, 0.97,
            (f"N = {n}\nr = {r:.2f}\nslope = {slope:.2f}\n"
             f"MB = {mb:.2e}\nRMSE = {rmse:.2e}"),
            transform=ax.transAxes, va="top", ha="right",
            fontsize=tk["fontsize"] * 0.8,
            bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.85))

    # make_scatter_density_plot's scatter branch parses "(units)" out of the
    # ylabel for its colorbar -- guarantee the parentheses exist.
    yl = ylabel or b["name"]
    units = yl[yl.find("(") + 1: yl.find(")")] if "(" in yl else "value"
    if "(" not in yl:
        yl = f"{b['name']} ({units})"
    xl = xlabel or (f"{a['name']} ({units})" if "(" not in (xlabel or a["name"])
                    else a["name"])

    df = pd.DataFrame({"x": xp, "y": yp})
    splots.make_scatter_density_plot(
        df, mod_var="x", obs_var="y", ax=ax, fill=bool(density),
        xlabel=xl, ylabel=yl, title=title,
        vmin_x=lo, vmax_x=hi, vmin_y=lo, vmax_y=hi,
        outname=outname + ".png")
    print(f"obs2obs: wrote {outname}.png (scatter, N={n})", flush=True)

def _series_points(spec, paired, tm, max_points=2_000_000):
    """Tidy per-record DataFrame [obs, mod, longitude, latitude, time,
    local_hour] for one obs_var/mod_var series -- works for both gridded
    (satellite obsgrid, flattened) and point (surface) pairs, with clip,
    bounds, and the local-hour window already applied and NaN rows dropped.

    The common currency for the platform-aware plot types (taylor, diurnal,
    xplatform_scatter): each needs obs and model side-by-side per record,
    not the pre-reduced ts/flat/grid of _extract.
    """
    obj = _pair_obj(paired, spec["label"])
    # value-mode series (single ``var``, as in the emissions boxplots) map
    # both obs and mod to that one field, so callers can use df["obs"].
    if "var" in spec and "mod_var" not in spec:
        ov = mv = spec["var"]
    else:
        ov, mv = spec["obs_var"], spec["mod_var"]
        
    if isinstance(obj, xr.Dataset):
        # build via numpy ravel 
        o = obj[ov]
        m = obj[mv].broadcast_like(o) if obj[mv].dims != o.dims else obj[mv]
        cols = {"obs": np.asarray(o.values).ravel(),
                "mod": np.asarray(m.values).ravel()}
        
        for c in ("longitude", "latitude"):
            if c in obj.coords or c in obj.variables:
                cols[c] = np.asarray(obj[c].broadcast_like(o).values).ravel()
        if "time" in o.coords or "time" in obj.coords:
            tt = obj["time"] if "time" in obj.coords else o["time"]
            cols["time"] = np.asarray(tt.broadcast_like(o).values).ravel()
        df = pd.DataFrame(cols)
    else:
        df = _as_dataframe(obj).rename(columns={ov: "obs", mv: "mod"})
        keep = ["obs", "mod"] + [c for c in ("longitude", "latitude", "time")
                                 if c in df.columns]
        df = df[keep].copy()

    # drop non-finite obs/mod BEFORE the rest -> the frame is small from here on
    df = df[np.isfinite(df["obs"]) & np.isfinite(df["mod"])]
    if "time" in df.columns:
        df["time"] = pd.to_datetime(df["time"])
        
    clip = _clip2(spec.get("clip"))
    if clip is not None:
        df = df[df["obs"].between(*clip) & df["mod"].between(*clip)]

    bnds = spec.get("bounds")
    if bnds is not None and "longitude" in df.columns:
        w, e, s, n = bnds
        df = df[(df.longitude >= w) & (df.longitude <= e)
                & (df.latitude >= s) & (df.latitude <= n)]

    if "longitude" in df.columns and "time" in df.columns:
        df["local_hour"] = ((df.time.dt.hour + df.time.dt.minute / 60.0
                             + df.longitude / 15.0) % 24.0)
        if isinstance(tm, dict) and "window_local" in tm:
            h0, h1 = tm["window_local"]
            df = df[(df.local_hour >= h0) & (df.local_hour <= h1)]

    df = df.dropna(subset=["obs", "mod"])
    if max_points and len(df) > max_points:
        df = df.iloc[:: int(np.ceil(len(df) / max_points))]
    return df.reset_index(drop=True)

def taylor(specs, paired, tm, outname, title=None, text_dict=None,
           fig_dict=None, max_points=2_000_000):
    """Normalized Taylor diagram: one point per platform (surface/TEMPO/
    TROPOMI). Each platform's std is normalized by its OWN obs std, so the
    reference std is 1 for all and different-unit platforms share one plot.
    """
    from monet.plots.taylordiagram import TaylorDiagram as td
    tk = {"fontsize": 14, **(text_dict or {})}
    _markers = ["o", "s", "^", "D", "v", "P", "*"]
    
    fig = plt.figure(figsize=(fig_dict or {}).get("figsize", (9, 8)))
    dia, drawn = None, 0
    for i, spec in enumerate(specs):
        try:
            df = _series_points(spec, paired, tm, max_points=max_points)
        except KeyError as e:
            print(f"obs2obs taylor: {e}; skipping.", flush=True)
            continue
        o, m = df["obs"].to_numpy(float), df["mod"].to_numpy(float)
        if o.size < 3 or o.std() == 0:
            print(f"obs2obs taylor: {spec.get('name', spec['label'])}: "
                  "too few points / zero variance; skipping.", flush=True)
            continue
        cc = float(np.corrcoef(o, m)[0, 1])
        nstd = float(m.std() / o.std())          # model std / obs std
        if dia is None:
            dia = td(1.0, scale=1.6, fig=fig, label="obs (ref)")
        dia.add_sample(nstd, cc, marker=_markers[i % len(_markers)], ms=11,
                       ls="", label=spec.get("name", spec["label"]),
                       color=spec.get("color"))
        drawn += 1
    if not drawn:
        plt.close(fig)
        return
    fig.legend(loc="upper right", fontsize=tk["fontsize"] * 0.8)
    if title:
        fig.suptitle(title, fontweight="bold", fontsize=tk["fontsize"])
    _save(outname + ".png")
    plt.close(fig)
    print(f"obs2obs: wrote {outname}.png (taylor, {drawn} platforms)", flush=True)

def diurnal(specs, paired, tm, outname, ylabel=None, title=None,
            normalize=False, text_dict=None, fig_dict=None,
            max_points=2_000_000):
    """Mean value vs local solar hour, obs and model, one line pair per
    series. ``normalize: true`` divides each curve by its own 24-h mean so
    different-unit platforms (surface ppb vs column) can share the axis as
    diurnal SHAPES. The window in ``tm`` is ignored here , so pass time_match: none/daily at the group.
    """
    tk = {"fontsize": 14, **(text_dict or {})}
    tm_full = "none" if isinstance(tm, dict) else tm    # never window a diurnal
    fig, ax = plt.subplots(figsize=(fig_dict or {}).get("figsize", (10, 6)))
    drawn = 0
    for spec in specs:
        try:
            df = _series_points(spec, paired, tm_full, max_points=max_points)
        except KeyError as e:
            print(f"obs2obs diurnal: {e}; skipping.", flush=True)
            continue
        if "local_hour" not in df or df.empty:
            continue
        hb = df["local_hour"].round().astype(int) % 24

        c = spec.get("color")
        nm = spec.get("name", spec["label"])
        if "var" in spec and "mod_var" not in spec:      # value mode: one line
            g = df.groupby(hb)["obs"].mean()
            if normalize and g.mean():
                g = g / g.mean()
            ax.plot(g.index, g.values, "-o", color=c, label=nm)
        else:                                            # obs + model lines
            go = df.groupby(hb)["obs"].mean()
            gm = df.groupby(hb)["mod"].mean()
            if normalize:
                go = go / go.mean() if go.mean() else go
                gm = gm / gm.mean() if gm.mean() else gm
            ax.plot(go.index, go.values, "-o", color=c, label=f"{nm} obs")
            ax.plot(gm.index, gm.values, "--s", color=c, alpha=0.7,
                    label=f"{nm} model")
        drawn += 1
    if not drawn:
        plt.close(fig)
        return
    ax.set_xlabel("local solar hour", fontweight="bold", fontsize=tk["fontsize"])
    ax.set_ylabel(ylabel or ("normalized (/24h mean)" if normalize else "value"),
                  fontweight="bold", fontsize=tk["fontsize"])
    ax.set_xticks(range(0, 24, 3))
    ax.legend(fontsize=tk["fontsize"] * 0.7, ncol=2)
    if title:
        ax.set_title(title, fontweight="bold", fontsize=tk["fontsize"])
    plt.tight_layout()
    _save(outname + ".png")
    plt.close(fig)
    print(f"obs2obs: wrote {outname}.png (diurnal, {drawn} series)", flush=True)
    
def xplatform_scatter(grp, paired, tm, outname, text_dict=None, fig_dict=None,
                      max_points=2_000_000):
    """Surface<->column coupling: x = surface conc, y = colocated satellite
    column, one cloud for OBS and one for MODEL, with best-fit slopes.

    Colocation is daily: surface sites are daily-averaged; the satellite
    column is daily-averaged per cell; each site-day is matched to the
    nearest column cell that same day (cKDTree). A model that preserves the
    vertical coupling (mixing height + profile) reproduces the observed
    surface->column slope; a flatter/steeper model slope does not.
    """
    from scipy.spatial import cKDTree
    scfg, ccfg = grp["surface"], grp["column"]
    for cfg in (scfg, ccfg):
        if grp.get("bounds") is not None and "bounds" not in cfg:
            cfg["bounds"] = grp["bounds"]
        if grp.get("clip") is not None and "clip" not in cfg:
            cfg["clip"] = grp["clip"]
    try:
        sdf = _series_points(scfg, paired, tm, max_points=max_points)
        cdf = _series_points(ccfg, paired, tm, max_points=max_points)
    except KeyError as e:
        print(f"obs2obs xplatform_scatter: {e}; skipping.", flush=True)
        return
    if sdf.empty or cdf.empty or "longitude" not in cdf:
        print("obs2obs xplatform_scatter: empty surface or column frame.", flush=True)
        return
    
    sdf = sdf.assign(day=sdf.time.dt.floor("D"))
    cdf = cdf.assign(day=cdf.time.dt.floor("D"))
    # daily mean per surface site (unique lon/lat) and per column cell
    s_day = (sdf.groupby(["day", "longitude", "latitude"])[["obs", "mod"]]
             .mean().reset_index())
    c_day = (cdf.groupby(["day", "longitude", "latitude"])[["obs", "mod"]]
             .mean().reset_index())

    rows = []
    for day, cg in c_day.groupby("day"):
        sg = s_day[s_day.day == day]
        if sg.empty or cg.empty:
            continue
        tree = cKDTree(cg[["longitude", "latitude"]].to_numpy())
        d, idx = tree.query(sg[["longitude", "latitude"]].to_numpy(), k=1)
        near = cg.iloc[idx]
        rows.append(pd.DataFrame({
            "sfc_obs": sg["obs"].to_numpy(), "sfc_mod": sg["mod"].to_numpy(),
            "col_obs": near["obs"].to_numpy(), "col_mod": near["mod"].to_numpy(),
            "sep_deg": d}))
    if not rows:
        print("obs2obs xplatform_scatter: no colocated site-days.", flush=True)
        return
    m = pd.concat(rows, ignore_index=True)
    m = m[m.sep_deg < grp.get("max_sep_deg", 0.5)]                # drop far matches
    m = m.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["sfc_obs", "sfc_mod", "col_obs", "col_mod"])
    if len(m) < 3:
        print("obs2obs xplatform_scatter: <3 colocated pairs.", flush=True)
        return
        
    tk = {"fontsize": 14, **(text_dict or {})}
    fig, ax = plt.subplots(figsize=(fig_dict or {}).get("figsize", (8, 7)))

    def _cloud(x, y, color, label):
        ax.scatter(x, y, s=10, alpha=0.4, color=color, edgecolors="none")
        sl, ic = np.polyfit(x, y, 1)
        r = np.corrcoef(x, y)[0, 1]
        xs = np.linspace(np.nanmin(x), np.nanmax(x), 50)
        ax.plot(xs, sl * xs + ic, color=color, lw=2,
                label=f"{label}: slope={sl:.2e}, r={r:.2f}")
        return sl, r

    _cloud(m.sfc_obs.to_numpy(), m.col_obs.to_numpy(), "black", "observed")
    _cloud(m.sfc_mod.to_numpy(), m.col_mod.to_numpy(), "purple", "model")
    ax.set_xlabel(grp.get("xlabel", f"surface {scfg['obs_var']}"),
                  fontweight="bold", fontsize=tk["fontsize"])
    ax.set_ylabel(grp.get("ylabel", f"column {ccfg['obs_var']}"),
                  fontweight="bold", fontsize=tk["fontsize"])
    ax.legend(fontsize=tk["fontsize"] * 0.75)
    ax.set_title(grp.get("title", "surface↔column coupling"),
                 fontweight="bold", fontsize=tk["fontsize"])
    plt.tight_layout()
    _save(outname + ".png")
    plt.close(fig)
    print(f"obs2obs: wrote {outname}.png (xplatform_scatter, N={len(m)})", flush=True)

def norm_spatial(grp, paired, tm, outname, text_dict=None, fig_dict=None):
    """Normalized multiscale map: satellite column time-mean field + surface
    site time-means, each divided by its domain mean (dimensionless ratio),
    obs and model side by side on a shared 0-2 scale. Shows whether surface
    and column hotspots co-locate and whether the model reproduces both.
    """
    scfg, ccfg = grp["surface"], grp["column"]
    cobj = _pair_obj(paired, ccfg["label"])
    sobj = _pair_obj(paired, scfg["label"])
    bnds = grp.get("bounds")
    tk = {"fontsize": 14, **(text_dict or {})}

    def _point_df(obj, var):
        """[longitude, latitude, var] rows, whether the surface pair is a
        pandas DataFrame or an xarray Dataset (netcdf-backed pairs load as
        the latter)."""
        if isinstance(obj, xr.Dataset):
            ds = obj[[var]]
            for c in ("longitude", "latitude"):
                if c in obj.coords or c in obj.variables:
                    ds = ds.assign_coords({c: obj[c]})
            d = ds.to_dataframe().reset_index()
        else:
            d = _as_dataframe(obj)
        return d[["longitude", "latitude", var]].dropna()
        
    def _field(var):
        da = cobj[var]
        da = da.mean("time", skipna=True) if "time" in da.dims else da
        mean = float(np.nanmean(da.values))
        return da / mean if mean else da

    def _sites(var):
        d = _point_df(sobj, var)
        if bnds is not None:
            w, e, s, n = bnds
            d = d[(d.longitude >= w) & (d.longitude <= e)
                  & (d.latitude >= s) & (d.latitude <= n)]
        g = d.groupby(["longitude", "latitude"])[var].mean().reset_index()
        mean = g[var].mean()
        g["r"] = g[var] / mean if mean else g[var]
        return g
    
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig, axes = plt.subplots(1, 2, figsize=(fig_dict or {}).get("figsize", (16, 6)),
                             subplot_kw={"projection": ccrs.PlateCarree()})
    for ax, (fvar, svar, ttl) in zip(axes, [
            (ccfg["obs_var"], scfg["obs_var"], "observed"),
            (ccfg["mod_var"], scfg["mod_var"], "model")]):
        fld, sit = _field(fvar), _sites(svar)
        pm = ax.pcolormesh(cobj["longitude"], cobj["latitude"], fld,
                           vmin=0, vmax=2, cmap="RdBu_r", shading="auto",
                           transform=ccrs.PlateCarree())
        ax.scatter(sit.longitude, sit.latitude, c=sit["r"], vmin=0, vmax=2,
                   cmap="RdBu_r", s=60, edgecolors="k", linewidths=0.6,
                   transform=ccrs.PlateCarree(), zorder=5)
        ax.coastlines(linewidth=0.5)
        ax.add_feature(cfeature.STATES, linewidth=0.3)
        if bnds is not None:
            ax.set_extent(bnds, crs=ccrs.PlateCarree())
        ax.set_title(ttl, fontweight="bold", fontsize=tk["fontsize"])
    cb = fig.colorbar(pm, ax=axes, orientation="horizontal", shrink=0.6,
                      pad=0.05, extend="both")
    cb.set_label("value / domain mean (column field, surface dots)",
                 fontweight="bold", fontsize=tk["fontsize"])
    if grp.get("title"):
        fig.suptitle(grp["title"], fontweight="bold", fontsize=tk["fontsize"])
    _save(outname + ".png")
    plt.close(fig)
    print(f"obs2obs: wrote {outname}.png (norm_spatial)", flush=True)


# run

_PLOTTERS = {"multi_boxplot": multi_boxplot,
             "multi_timeseries": multi_timeseries,
             "diff_map": diff_map,
             "scatter": scatter}

# plot types that take the whole group because they need their own obs/model alignment
_GROUP_PLOTTERS = {"taylor": taylor, "diurnal": diurnal,
                   "xplatform_scatter": xplatform_scatter,
                   "norm_spatial": norm_spatial}

# spatial overlay 
# spatial bias 
# 

def labels_used(config, only=None):
    """Every pair label referenced by the obs2obs groups (optionally only the
    groups whose name contains ``only``). Use this to prune
    read.paired.filenames BEFORE read_analysis so a run opens only the files
    it needs 
    """
    used = set()
    for gname, grp in config.get("groups", {}).items():
        if only and only not in gname:
            continue
        for s in grp.get("series", []):
            if s.get("label"):
                used.add(s["label"])
        for key in ("surface", "column"):        # xplatform/norm group configs
            if isinstance(grp.get(key), dict) and grp[key].get("label"):
                used.add(grp[key]["label"])
    return used
    
def run(paired, config, default_outdir=".", only = None):
    """Execute every group in the ``obs2obs:`` config against loaded pairs."""
    outdir = os.path.expandvars(config.get("output_dir", default_outdir))
    os.makedirs(outdir, exist_ok=True)
    default_tm = config.get("time_match", "none")
    max_points = int(config.get("max_points", 2_000_000))

    for gname, grp in config.get("groups", {}).items():
        if only and only not in gname:
            continue
        dp = grp.get("data_proc") or {}
        tm = grp.get("time_match", default_tm)
        gbounds = grp.get("bounds")
        gclip = grp.get("clip")

        # group-level plot types (taylor/diurnal/xplatform_scatter/norm_spatial)
        # handle their own extraction/colocation from the raw group config.
        types = grp.get("type", "multi_boxplot")
        types = [types] if isinstance(types, str) else list(types)
        _group_types = [t for t in types if t in _GROUP_PLOTTERS]
        
        for t in _group_types:
            out = os.path.join(outdir, f"obs2obs_{gname}.{t}")
            kw = dict(text_dict=grp.get("text_kwargs"),
                      fig_dict=grp.get("fig_kwargs"), max_points=max_points)
            try:
                if t == "taylor":
                    _GROUP_PLOTTERS[t](grp["series"], paired, tm, out,
                                       title=grp.get("title", gname), **kw)
                elif t == "diurnal":
                    _GROUP_PLOTTERS[t](grp["series"], paired, tm, out,
                                       ylabel=grp.get("ylabel"),
                                       title=grp.get("title", gname),
                                       normalize=dp.get("normalize",
                                                         grp.get("normalize", False)),
                                       **kw)
                else:                       # xplatform_scatter / norm_spatial
                    _kw = {k: v for k, v in kw.items() if k != "max_points"}
                    if t == "xplatform_scatter":
                        _kw["max_points"] = max_points
                    _GROUP_PLOTTERS[t](grp, paired, tm, out, **_kw)
            except Exception as e:  # noqa: BLE001
                print(f"obs2obs: {gname}.{t} failed: {e!r}", flush=True)
        types = [t for t in types if t not in _GROUP_PLOTTERS]
        if not types:
            continue
            
        series = []
        for spec in grp["series"]:
            if gbounds is not None and "bounds" not in spec:
                spec = {**spec, "bounds": gbounds}
            if gclip is not None and "clip" not in spec:
                spec = {**spec, "clip": gclip}
            try:
                series.append(_extract(spec, paired, tm, max_points=max_points))
            except KeyError as e:
                print(f"obs2obs: {gname}: {e}; skipping series.", flush=True)
        if not series:
            continue
        types = grp.get("type", "multi_boxplot")
        for t in types:
            fn = _PLOTTERS.get(t)
            if fn is None:
                print(f"obs2obs: unknown type {t!r} "
                      f"(use {sorted(_PLOTTERS)}).", flush=True)
                continue
            out = os.path.join(outdir, f"obs2obs_{gname}.{t}")
            kwargs = dict(ylabel=grp.get("ylabel"),
                          title=grp.get("title", gname),
                          fig_dict=grp.get("fig_kwargs"),
                          text_dict=grp.get("text_kwargs"))
            if t == "multi_boxplot":
                kwargs.update(
                    vmin=_fnum(grp.get("vmin")), vmax=_fnum(grp.get("vmax")),
                    showfliers=dp.get("showfliers",
                                      grp.get("showfliers", False)),
                    hline=_fnum(dp.get("hline", grp.get("hline"))))
            if t == "multi_timeseries":
                kwargs["avg_window"] = grp.get("avg_window")
            if t == "diff_map":
                kwargs.update(vdiff=_fnum(grp.get("vdiff")),
                              domain_name=grp.get("domain_name"))
            if t == "scatter":
                kwargs.update(xlabel=grp.get("xlabel"),
                              xylim=grp.get("xylim"),
                              density=dp.get("density", grp.get("density", False)),
                              time_avg=grp.get("time_avg", "1D"))
                
            fn(series, out, **kwargs)
        # free this group's extracted arrays before moving on
        del series
        gc.collect()

