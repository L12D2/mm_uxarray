
"""Native-resolution TROPOMI city pairing -- TROPOMI counterpart of
pair_tempo_native.py.

Env (all optional except RUN):
    RUN=nonbiog|biog|grapes|mxcat        (required)
    CITY=atl|mex|la|den|dfw              (default: all cities for the run)
    YMD=20240601                         (required; one day, PBS array sets it)
    OBS_GRID_RES=0.05                    (default 0.05 deg; TROPOMI ~5.5x3.5 km)
    REGRID_METHOD=conservative           (default conservative)
    REGRID_TARGET=obs|model|swath        (default obs)
    TROPOMI_PRODUCTS=no2,hcho,co         (default: all three)

Submit (mirrors submit_pair_tempo_native.sh -- set YMD from PBS_ARRAY_INDEX):
    for RUN in nonbiog biog grapes; do
      for CITY in atl dfw la den; do
        qsub -N tn_${RUN}_${CITY} -o tn_${RUN}_${CITY}.log \
             -l select=1:ncpus=1:mem=48GB \
             -v RUN=$RUN,CITY=$CITY,REGRID_METHOD=conservative,REGRID_TARGET=swath \
             submit_pair_tropomi_native.sh
      done
    done
"""

import os, glob, copy, time
from datetime import datetime, timedelta

import yaml
from melodies_monet import driver

# reuse the model definitions, CONUS extent, and method tags from the TEMPO driver
from pair_tempo_native import BASE, CONUS_EXTENT, METHOD_TAG, RUNS
from mm_paths import (paired_dir, gridtype_of, apply_filters, apply_save,
                      filter_tag, save_suffix, _envflag)

# TROPOMI L2 granule sources (one file per overpass; the S5P filename carries the
# start datetime YYYYMMDDTHHMMSS, so a per-day glob is "..._{ymd}T*"). Note the
# product field widths differ: NO2____ (4), HCHO___ (3), CO_____ (5).
OBS_SOURCES = {
    "tropomi_l2_no2":  {"dir": "/glade/derecho/scratch/lcthompson/tropomi/no2_2024",
                        "glob": "S5P_OFFL_L2__NO2____{ymd}T*.nc"},
    "tropomi_l2_hcho": {"dir": "/glade/derecho/scratch/lcthompson/tropomi/hcho_2024",
                        "glob": "S5P_OFFL_L2__HCHO___{ymd}T*.nc"},
    "tropomi_l2_co":   {"dir": "/glade/derecho/scratch/lcthompson/tropomi/co_2024",
                        "glob": "S5P_OFFL_L2__CO_____{ymd}T*.nc"},
}

# model species (mapping) per product -- from control_master.yaml model block
MAPPING = {
    "tropomi_l2_no2":  {"NO2":  "nitrogendioxide_tropospheric_column"},
    "tropomi_l2_hcho": {"CH2O": "formaldehyde_tropospheric_vertical_column"},
    "tropomi_l2_co":   {"CO":   "carbonmonoxide_total_column"},
}

# obs `variables` blocks -- what the reader must load (column + AK + geometry +
# per-pixel precision for uncertainty). QA screening is injected by
# mm_paths.apply_filters (qa_value always on; cloud/SZA per env toggles).
OBS_BLOCKS = {
    "tropomi_l2_no2": {
        "obs_type": "sat_swath_clm", "sat_type": "tropomi_l2_no2",
        "obs_grid_units": "deg",
        "variables": {
            "qa_value": {},
            "nitrogendioxide_tropospheric_column": {},
            "averaging_kernel": {},
            "air_mass_factor_troposphere": {},
            "air_mass_factor_total": {},
            "tm5_tropopause_pressure": {},
            "latitude_bounds": {}, "longitude_bounds": {},
            "nitrogendioxide_tropospheric_column_precision": {},
        },
    },
    "tropomi_l2_hcho": {
        "obs_type": "sat_swath_clm", "sat_type": "tropomi_l2_hcho",
        "obs_grid_units": "deg",
        "variables": {
            "qa_value": {},
            "formaldehyde_tropospheric_vertical_column": {},
            "averaging_kernel": {},
            "formaldehyde_tropospheric_air_mass_factor": {},
            "surface_pressure": {},
            "tm5_constant_a": {}, "tm5_constant_b": {},
            "latitude_bounds": {}, "longitude_bounds": {},
            "formaldehyde_tropospheric_vertical_column_precision": {},
        },
    },
    "tropomi_l2_co": {
        "obs_type": "sat_swath_clm", "sat_type": "tropomi_l2_co",
        "obs_grid_units": "deg",
        "variables": {
            "qa_value": {},
            "carbonmonoxide_total_column": {},
            "column_averaging_kernel": {},
            "pressure_levels": {},
            "latitude_bounds": {}, "longitude_bounds": {},
            "carbonmonoxide_total_column_precision": {},
            "carbonmonoxide_total_column_corrected": {},
        },
    },
}

_SHORT = {"no2": "tropomi_l2_no2", "hcho": "tropomi_l2_hcho", "co": "tropomi_l2_co"}


def main():
    ymd = os.environ["YMD"]
    run = os.environ["RUN"].strip()
    if run not in RUNS:
        raise SystemExit(f"RUN={run!r} not in {list(RUNS)}")
    res = float(os.environ.get("OBS_GRID_RES", "0.05"))
    restag = ("%g" % res).replace(".", "p")
    method = os.environ.get("REGRID_METHOD", "conservative").strip()
    methodtag = METHOD_TAG.get(method, method[:5])
    target = os.environ.get("REGRID_TARGET", "model").strip().lower()
    if target not in ("obs", "model", "swath"):
        raise SystemExit(f"REGRID_TARGET={target!r} must be 'obs', 'model', or 'swath'")
    ftag = filter_tag(_envflag("CLOUD_FILTER"), _envflag("SZA_FILTER")) + save_suffix()

    want = os.environ.get("TROPOMI_PRODUCTS", "no2,hcho,co")
    products = [_SHORT[s.strip()] for s in want.split(",") if s.strip() in _SHORT]
    if not products:
        raise SystemExit(f"TROPOMI_PRODUCTS={want!r} resolved to nothing")

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

    # mm_output_v2/<run-slug>/tropomi/<gridtype>  (model | swath | grid_0p05 ...)
    gridtype = gridtype_of(target, restag)
    outdir = paired_dir(run, "tropomi", gridtype)
    yamldir = f"{outdir}/yaml"
    os.makedirs(yamldir, exist_ok=True)

    # HOURS narrows the granule glob for cheap smoke tests (all glob templates
    # contain "{ymd}T"; splice the hour pattern right after that T).
    _hh = os.environ.get("HOURS", "").strip()
    obs_pat = {}
    for name in products:
        _p = f"{OBS_SOURCES[name]['dir']}/{OBS_SOURCES[name]['glob'].format(ymd=ymd)}"
        obs_pat[name] = _p.replace(f"{ymd}T", f"{ymd}T{_hh}") if _hh else _p
    if _hh:
        print(f"HOURS={_hh!r}: granule glob narrowed to {ymd}T{_hh}*", flush=True)
    have = {name: bool(glob.glob(pat)) for name, pat in obs_pat.items()}
    if not any(have.values()):
        print(f"SKIP {run} {ymd}: no TROPOMI granules ({list(obs_pat.values())})",
              flush=True)
        return

    t0 = time.time()

    cd = copy.deepcopy(yaml.safe_load(open(BASE)))
    cd["analysis"]["start_time"] = iso
    cd["analysis"]["end_time"] = f"{iso} 23:59:00"
    for k in ("output_dir", "output_dir_save", "output_dir_read"):
        cd["analysis"][k] = outdir
    cd["model"]["cam-chem-se"]["files"] = mfiles
    cd["model"]["cam-chem-se"]["scrip_file"] = r["scrip"]

    cd["obs"] = {}
    mapping = {}
    present = []
    for name in products:
        if not have.get(name):
            continue
        blk = copy.deepcopy(OBS_BLOCKS[name])
        blk["filename"] = obs_pat[name]
        blk["regrid_method"] = method
        blk["regrid_target"] = [target]
        # ALWAYS set the extent: it is the granule crop for every target
        # (a stale extent in the BASE control would silently crop the run).
        blk["obs_grid_extent"] = CONUS_EXTENT
        if target == "obs":                # res only builds the lat/lon grid
            blk["obs_grid_res"] = res
        apply_filters(blk, name)           # QA screen per CLOUD/SZA env
        apply_save(blk, name) 
        cd["obs"][name] = blk
        mapping[name] = MAPPING[name]
        present.append(name)
    cd["model"]["cam-chem-se"]["mapping"] = mapping
    if not present:
        print(f"SKIP {run} {ymd}: no obs products with data", flush=True)
        return

    _grid_bit = f"conus_{restag}_" if target == "obs" else f"{target}_"
    prefix = f"{_grid_bit}{methodtag}_{ftag}_{ymd}"
    cd["analysis"]["save"]["paired"]["prefix"] = prefix

    tmp = f"{yamldir}/control_{prefix}.yaml"
    yaml.safe_dump(cd, open(tmp, "w"), sort_keys=False)
    print(f"==== {run} {iso} CONUS [{method}->{target}] {ftag} "
          f"products={present} ====", flush=True)
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
