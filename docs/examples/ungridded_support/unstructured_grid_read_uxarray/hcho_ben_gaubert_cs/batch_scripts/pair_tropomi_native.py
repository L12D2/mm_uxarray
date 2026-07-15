
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

# reuse the model definitions, city boxes, and method tags from the TEMPO driver
from pair_tempo_native import BASE, CITIES, METHOD_TAG, RUNS

OUTROOT = "/glade/work/lcthompson/mm_output"

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
# per-pixel precision for uncertainty). Lifted from the user's control YAML.
OBS_BLOCKS = {
    "tropomi_l2_no2": {
        "obs_type": "sat_swath_clm", "sat_type": "tropomi_l2_no2",
        "obs_grid_units": "deg",
        "variables": {
            "qa_value": {"qa_min": 0.75},
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
            "qa_value": {"qa_min": 0.5},
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
            "qa_value": {"qa_min": 0.5},
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
    target = os.environ.get("REGRID_TARGET", "obs").strip().lower()
    if target not in ("obs", "model", "swath"):
        raise SystemExit(f"REGRID_TARGET={target!r} must be 'obs', 'model', or 'swath'")

    want = os.environ.get("TROPOMI_PRODUCTS", "no2,hcho,co")
    products = [_SHORT[s.strip()] for s in want.split(",") if s.strip() in _SHORT]
    if not products:
        raise SystemExit(f"TROPOMI_PRODUCTS={want!r} resolved to nothing")

    city_env = os.environ.get("CITY", "").strip().lower()
    if city_env:
        if city_env not in CITIES:
            raise SystemExit(f"CITY={city_env!r} not in {list(CITIES)}")
        cities = {city_env: CITIES[city_env]}
    else:
        cities = CITIES

    iso = datetime.strptime(ymd, "%Y%m%d").strftime("%Y-%m-%d")
    r = RUNS[run]
    if not glob.glob(f"{r['model_dir']}/{r['model_stem']}.{iso}-*.nc"):
        print(f"SKIP {run} {ymd}: no model files under {r['model_dir']}", flush=True)
        return

    days = [(d0 + timedelta(days=k)).strftime("%Y-%m-%d") for k in (-1, 0, 1)]
    mfiles = sorted(f for day in days
                    for f in glob.glob(f"{r['model_dir']}/{r['model_stem']}.{day}-*.nc"))
    
    outdir = f"{OUTROOT}/{run}_tropomi_native"
    yamldir = f"{outdir}/yaml"
    os.makedirs(yamldir, exist_ok=True)

    obs_pat = {name: f"{OBS_SOURCES[name]['dir']}/"
                     f"{OBS_SOURCES[name]['glob'].format(ymd=ymd)}"
               for name in products}
    have = {name: bool(glob.glob(pat)) for name, pat in obs_pat.items()}
    if not any(have.values()):
        print(f"SKIP {run} {ymd}: no TROPOMI granules ({list(obs_pat.values())})",
              flush=True)
        return

    t0 = time.time()

    def _base_cd():
        cd = copy.deepcopy(yaml.safe_load(open(BASE)))
        cd["analysis"]["start_time"] = iso
        cd["analysis"]["end_time"] = f"{iso} 23:59:00"
        for k in ("output_dir", "output_dir_save", "output_dir_read"):
            cd["analysis"][k] = outdir
        cd["model"]["cam-chem-se"]["files"] = mfiles
        cd["model"]["cam-chem-se"]["scrip_file"] = r["scrip"]
        return cd

    def _setup_obs(cd, box=None):
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
            if target in ("obs", "swath"):     # 'model' ignores res/extent
                blk["obs_grid_res"] = res
                blk["obs_grid_extent"] = box
            cd["obs"][name] = blk
            mapping[name] = MAPPING[name]
            present.append(name)
        cd["model"]["cam-chem-se"]["mapping"] = mapping
        return present

    def _run(cd, tmp, tag):
        yaml.safe_dump(cd, open(tmp, "w"), sort_keys=False)
        print(f"==== {run} {iso} {tag} [{method}->{target}] ====", flush=True)
        an = driver.analysis()
        an.control = tmp
        an.read_control()
        an.open_models()
        an.open_obs()
        an.pair_data()
        an.save_analysis()
        for label, p in an.paired.items():
            print(f"  {tag}: {label}: {dict(p.obj.sizes)}", flush=True)

    if target == "model":
        # full-domain mesh product; city extent/res do not apply
        cd = _base_cd()
        cd["analysis"]["save"]["paired"]["prefix"] = f"model_{methodtag}_{ymd}"
        if _setup_obs(cd):
            _run(cd, f"{yamldir}/control_model_{methodtag}_{ymd}.yaml", "model-space")
    else:
        for city, box in cities.items():
            cd = _base_cd()
            cd["analysis"]["save"]["paired"]["prefix"] = f"{city}_{restag}_{methodtag}_{ymd}"
            if _setup_obs(cd, box):
                _run(cd, f"{yamldir}/control_{city}_{restag}_{methodtag}_{ymd}.yaml",
                     f"{city}@{res}deg")

    print(f"==== DONE {run} {iso} in {time.time()-t0:.0f}s ====", flush=True)


if __name__ == "__main__":
    main()