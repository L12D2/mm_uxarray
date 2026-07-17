"""CONUS-scale TEMPO pairing, one day per PBS array element (env YMD).

Env (all optional except RUN, YMD):
    RUN=nonbiog|biog|grapes|mxcat       (required)
    YMD=20240608                        (required; set from PBS_ARRAY_INDEX)
    REGRID_TARGET=model|obs|swath       (default model; obs = CONUS lat/lon grid)
    REGRID_METHOD=conservative|radius_mean|...   (default conservative)
    OBS_GRID_RES=0.03                   (obs target only; deg)
    CLOUD_FILTER=on  SZA_FILTER=on      (QA screen toggles; base flag always on)
    HOURS='1[5-7]'                      (smoke test: only {ymd}T15*..T17* granules)

Everything is CONUS-wide: granules are cropped to CONUS_EXTENT for every
target. City-level zooms were removed -- zoom at PLOT time instead.
"""

import os, glob, copy, time
from datetime import datetime, timedelta
import yaml
from melodies_monet import driver
from mm_paths import (paired_dir, gridtype_of, apply_filters, apply_save,
                      filter_tag, save_suffix, _envflag)

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "control_tempo_native.yaml")

_NE0 = "/glade/campaign/acom/MUSICA/grids/ne0CONUSne30x8/ne0CONUS_ne30x8_np4_SCRIP.nc"

# model month source + grid per emissions run (H1 hourly, for overpass matching)
RUNS = {
    "nonbiog": dict(
        model_dir="/glade/campaign/acom/acom-da/conus_outputs/"
                  "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.ERA5_ref_dust_M1.1.002/H1",
        model_stem="f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.ERA5_ref_dust_M1.1.002.cam.h1",
        scrip=_NE0),
    "biog": dict(
        model_dir="/glade/campaign/acom/acom-da/conus_outputs/"
                  "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.biog_ERA5_ref_dust_M1.1.001/H1",
        model_stem="f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.biog_ERA5_ref_dust_M1.1.001.cam.h1",
        scrip=_NE0),
    "grapes": dict(
        model_dir="/glade/campaign/acom/acom-da/conus_outputs/"
                  "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.nox_grapes.001/H1",
        model_stem="f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.nox_grapes.001.cam.h1",
        scrip=_NE0),
    "mxcat": dict(
        model_dir="/glade/campaign/acom/acom-weather/jjacdan/SCENICS.HAMAQ/"
                  "f.e3beta01.FCts2nudged.MXCATL_ne30x16_cams_mosaic_v1.1_final.03",
        model_stem="f.e3beta01.FCts2nudged.MXCATL_ne30x16_cams_mosaic_v1.1_final.03.cam.h3i",
        scrip="/glade/work/jjacdan/ne0np4.MXC.ATL.ne30x16/grids/MXC.ATL_ne30x16_np4_SCRIP.nc"),
}

# City boxes: retained ONLY because pair_tropomi_native imports CITIES; the
# city pairing path was removed from this script (CONUS-scale only).
CITIES = {
    "atl":  [-85.0, -83.7, 33.2, 34.4],
    "mex":  [-99.6, -98.6, 18.9, 20.0],
    "la":   [-119.2, -117.3, 33.3, 34.7],
    "den":  [-105.8, -104.3, 39.2, 40.3],
    "dfw":  [-97.8, -96.3, 32.2, 33.4],
}

# CONUS-wide crop [lonW, lonE, latS, latN] (latS=15 keeps Mexico City / TEMPO's
# southern FOR edge)
CONUS_EXTENT = [-130.0, -60.0, 15.0, 55.0]

OBS_SOURCES = {
    "tempo_l2_hcho": {"dir": "/glade/campaign/acom/acom-da/sma/TEMPO_HCHO_V03",
                      "glob": "TEMPO_HCHO_L2_V03_{ymd}T*_S*"},
    "tempo_l2_no2":  {"dir": "/glade/campaign/acom/acom-da/sma/TEMPO_NO2_V03",
                      "glob": "TEMPO_NO2_L2_V03_{ymd}T*_S*"},
}

METHOD_TAG = {"conservative": "cons", "conservative_normed": "consn",
              "radius_mean": "rmean", "nearest_s2d": "nn", "nearest_d2s": "nnd",
              "bilinear": "bilin"}


def main():
    ymd = os.environ["YMD"]
    run = os.environ["RUN"].strip()
    if run not in RUNS:
        raise SystemExit(f"RUN={run!r} not in {list(RUNS)}")
    res = float(os.environ.get("OBS_GRID_RES", "0.03"))
    restag = ("%g" % res).replace(".", "p")
    method = os.environ.get("REGRID_METHOD", "conservative").strip()
    methodtag = METHOD_TAG.get(method, method[:5])
    target = os.environ.get("REGRID_TARGET", "model").strip().lower()
    if target not in ("obs", "model", "swath"):
        raise SystemExit(f"REGRID_TARGET={target!r} must be 'obs', 'model', or 'swath'")
    # QA-screen sensitivity toggle (CLOUD_FILTER / SZA_FILTER env); base flag always on
    ftag = filter_tag(_envflag("CLOUD_FILTER"), _envflag("SZA_FILTER")) + save_suffix()

    # Running
    # one product per job halves the per-day granule accumulation 
    _want = os.environ.get("TEMPO_PRODUCTS", "hcho,no2")
    keep = {f"tempo_l2_{s.strip()}" for s in _want.split(",")
            if s.strip() in ("hcho", "no2")}
    if not keep:
        raise SystemExit(f"TEMPO_PRODUCTS={_want!r} resolved to nothing")

    d0 = datetime.strptime(ymd, "%Y%m%d")
    iso = d0.strftime("%Y-%m-%d")
    r = RUNS[run]
    if not glob.glob(f"{r['model_dir']}/{r['model_stem']}.{iso}-*.nc"):
        print(f"SKIP {run} {ymd}: no model files under {r['model_dir']}", flush=True)
        return

    # neighbor days so 00 UTC is covered (CAM h1 end-of-interval stamping)
    days = [(d0 + timedelta(days=k)).strftime("%Y-%m-%d") for k in (-1, 0, 1)]
    mfiles = sorted(f for day in days
                    for f in glob.glob(f"{r['model_dir']}/{r['model_stem']}.{day}-*.nc"))

    # mm_output_v2/<run-slug>/tempo/<gridtype>  (model | swath | grid_0p03 ...)
    gridtype = gridtype_of(target, restag)
    outdir = paired_dir(run, "tempo", gridtype)
    yamldir = f"{outdir}/yaml"
    os.makedirs(yamldir, exist_ok=True)

    # obs granule patterns for this day. 
    _hh = os.environ.get("HOURS", "").strip()
    obs_pat = {name: f"{src['dir']}/{src['glob'].format(ymd=ymd)}".replace(
                   f"{ymd}T", f"{ymd}T{_hh}") if _hh
               else f"{src['dir']}/{src['glob'].format(ymd=ymd)}"
               for name, src in OBS_SOURCES.items() if name in keep}
    if _hh:
        print(f"HOURS={_hh!r}: granule glob narrowed to {ymd}T{_hh}*", flush=True)
    have = {name: bool(glob.glob(pat)) for name, pat in obs_pat.items()}
    if not any(have.values()):
        print(f"SKIP {run} {ymd}: no TEMPO granules", flush=True)
        return

    t0 = time.time()

    cd = copy.deepcopy(yaml.safe_load(open(BASE)))
    cd["analysis"]["start_time"] = iso
    cd["analysis"]["end_time"] = f"{iso} 23:59:00"
    for k in ("output_dir", "output_dir_save", "output_dir_read"):
        cd["analysis"][k] = outdir
    cd["model"]["cam-chem-se"]["files"] = mfiles
    cd["model"]["cam-chem-se"]["scrip_file"] = r["scrip"]

    mapping = cd["model"]["cam-chem-se"].get("mapping", {})
    present = []
    for name in list(OBS_SOURCES):
        if have.get(name):
            cd["obs"][name]["filename"] = obs_pat[name]
            cd["obs"][name]["regrid_target"] = [target]
            cd["obs"][name]["regrid_method"] = method
            # ALWAYS set the extent: it is the granule crop for every target.
            cd["obs"][name]["obs_grid_extent"] = CONUS_EXTENT
            if target == "obs":            # res only builds the lat/lon grid
                cd["obs"][name]["obs_grid_res"] = res
            apply_filters(cd["obs"][name], name)   # QA screen per CLOUD/SZA env
            apply_save(cd["obs"][name], name)
            present.append(name)
        else:
            cd["obs"].pop(name, None)
            mapping.pop(name, None)
    if not present:
        print(f"SKIP {run} {ymd}: no obs products with data", flush=True)
        return

    _grid_bit = f"conus_{restag}_" if target == "obs" else f"{target}_"
    prefix = f"{_grid_bit}{methodtag}_{ftag}_{ymd}"
    cd["analysis"]["save"]["paired"]["prefix"] = prefix

    tmp = f"{yamldir}/control_{prefix}.yaml"
    yaml.safe_dump(cd, open(tmp, "w"), sort_keys=False)
    print(f"==== {run} {iso} CONUS [{method}->{target}] {ftag} ====", flush=True)
    an = driver.analysis()
    an.control = tmp
    an.read_control()
    an.open_models()
    an.open_obs()
    an.pair_data()
    an.save_analysis()
    for label, p in an.paired.items():
        print(f"  {label}: {dict(p.obj.sizes)}", flush=True)

    print(f"==== DONE {run} {iso} in {time.time()-t0:.0f}s ====", flush=True)


if __name__ == "__main__":
    main()
