
import os, glob, copy, time
from datetime import datetime
import yaml
from melodies_monet import driver

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "control_tempo_native.yaml")
OUTROOT = "/glade/work/lcthompson/mm_output"

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

# [lonW, lonE, latS, latN]
CITIES = {
    "atl":  [-85.0, -83.7, 33.2, 34.4],
    "mex":  [-99.6, -98.6, 18.9, 20.0],
    #"seus": [-92.0, -75.0, 24.5, 37.0],   
    "la":   [-119.2, -117.3, 33.3, 34.7],
    "den":  [-105.8, -104.3, 39.2, 40.3],
    "dfw":  [-97.8, -96.3, 32.2, 33.4],
}

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
    target = os.environ.get("REGRID_TARGET", "obs").strip().lower()
    if target not in ("obs", "model"):
        raise SystemExit(f"REGRID_TARGET={target!r} must be 'obs' or 'model'")

    city_env = os.environ.get("CITY", "").strip().lower()
    if city_env:
        if city_env not in CITIES:
            raise SystemExit(f"CITY={city_env!r} not in {list(CITIES)}")
        cities = {city_env: CITIES[city_env]}
    else:
        cities = CITIES

    iso = datetime.strptime(ymd, "%Y%m%d").strftime("%Y-%m-%d")
    r = RUNS[run]
    mfiles = sorted(glob.glob(f"{r['model_dir']}/{r['model_stem']}.{iso}-*.nc"))
    if not mfiles:
        print(f"SKIP {run} {ymd}: no model files under {r['model_dir']}", flush=True)
        return

    outdir = f"{OUTROOT}/{run}_tempo_native"
    yamldir = f"{outdir}/yaml"
    os.makedirs(yamldir, exist_ok=True)

    # obs granule patterns for this day, resolved once
    obs_pat = {name: f"{src['dir']}/{src['glob'].format(ymd=ymd)}"
               for name, src in OBS_SOURCES.items()}
    have = {name: bool(glob.glob(pat)) for name, pat in obs_pat.items()}
    if not any(have.values()):
        print(f"SKIP {run} {ymd}: no TEMPO granules", flush=True)
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
        mapping = cd["model"]["cam-chem-se"].get("mapping", {})
        present = []
        for name in list(OBS_SOURCES):
            if have.get(name):
                cd["obs"][name]["filename"] = obs_pat[name]
                cd["obs"][name]["regrid_target"] = [target]
                cd["obs"][name]["regrid_method"] = method
                if target == "obs":                    # mesh target ignores these
                    cd["obs"][name]["obs_grid_res"] = res
                    cd["obs"][name]["obs_grid_extent"] = box
                present.append(name)
            else:
                cd["obs"].pop(name, None)
                mapping.pop(name, None)
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