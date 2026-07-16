"""
pair ONE day of TEMPO L2 HCHO vs CESM-SE. Driven by env YMD
"""

import os, glob, copy, time
from datetime import datetime, timedelta
import yaml
from melodies_monet import driver

BASE = ("/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/"
        "unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/control_tempo_l2_hcho_cesm_se.yaml")

OUTROOT = "/glade/work/lcthompson/mm_output"
_NE0 = "/glade/campaign/acom/MUSICA/grids/ne0CONUSne30x8/ne0CONUS_ne30x8_np4_SCRIP.nc"

# per-emissions-run model source + mesh + paired dir (mirrors pair_tempo_native.py
# and runs.yaml)
RUNS = {
    "nonbiog": dict(
        model_dir="/glade/campaign/acom/acom-da/conus_outputs/"
                  "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.ERA5_ref_dust_M1.1.002/H1",
        model_stem="f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.ERA5_ref_dust_M1.1.002.cam.h1",
        scrip=_NE0,
        paired_dir=f"{OUTROOT}/nonbiog_refera5_dust"),
    "biog": dict(
        model_dir="/glade/campaign/acom/acom-da/conus_outputs/"
                  "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.biog_ERA5_ref_dust_M1.1.001/H1",
        model_stem="f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.biog_ERA5_ref_dust_M1.1.001.cam.h1",
        scrip=_NE0,
        paired_dir=f"{OUTROOT}/biog_refera5_dust"),
    "grapes": dict(
        model_dir="/glade/campaign/acom/acom-da/conus_outputs/"
                  "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.nox_grapes.001/H1",
        model_stem="f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.nox_grapes.001.cam.h1",
        scrip=_NE0,
        paired_dir=f"{OUTROOT}/grapes"),
    "mxcat": dict(
        model_dir="/glade/campaign/acom/acom-weather/jjacdan/SCENICS.HAMAQ/"
                  "f.e3beta01.FCts2nudged.MXCATL_ne30x16_cams_mosaic_v1.1_final.03",
        model_stem="f.e3beta01.FCts2nudged.MXCATL_ne30x16_cams_mosaic_v1.1_final.03.cam.h3i",
        scrip="/glade/work/jjacdan/ne0np4.MXC.ATL.ne30x16/grids/MXC.ATL_ne30x16_np4_SCRIP.nc",
        paired_dir=f"{OUTROOT}/mxcat"),
}

RUN = os.environ.get("RUN", "nonbiog").strip()
if RUN not in RUNS:
    raise SystemExit(f"RUN={RUN!r} not in {sorted(RUNS)}")
_R = RUNS[RUN]

OUTDIR = _R["paired_dir"]
YAMLDIR = OUTDIR + "/hcho_yaml"
MODEL_DIR = _R["model_dir"]
MODEL_STEM = _R["model_stem"]

OBS_SOURCES = {
    
    "tempo_l2_hcho": {"dir": "/glade/campaign/acom/acom-da/sma/TEMPO_HCHO_V03",
                      "glob": "TEMPO_HCHO_L2_V03_{ymd}T*_S*"},
    "tempo_l2_no2": {"dir": "/glade/campaign/acom/acom-da/sma/TEMPO_NO2_V03",
                 "glob": "TEMPO_NO2_L2_V03_{ymd}T*_S*"},
    
    "tropomi_l2_hcho": {"dir": "/glade/derecho/scratch/lcthompson/tropomi/hcho_2024",
                        "glob": "S5P_*_L2__HCHO___{ymd}T*.nc"},
    "tropomi_l2_no2": {"dir": "/glade/derecho/scratch/lcthompson/tropomi/no2_2024",
                        "glob": "S5P_*_L2__NO2____{ymd}T*.nc"},

    "tropomi_l2_co": {"dir": "/glade/derecho/scratch/lcthompson/tropomi/co_2024",
                        "glob": "S5P_*_L2__CO_____{ymd}T*.nc"},
    
    # "tropomi_l2_hcho": {"dir": "/glade/campaign/acom/acom-da/sma/TROPOMI-HCHO-DATA/2024",
    #                 "glob": "S5P_OFFL_L2__HCHO___{ymd}T*.nc"},
}

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(YAMLDIR, exist_ok=True)

def main():
    ymd = os.environ["YMD"]
    d0 = datetime.strptime(ymd, "%Y%m%d")
    iso = d0.strftime("%Y-%m-%d")
    t0 = time.time()

    # Require the target day's model file (skip logic keys off THIS day only).
    if not glob.glob(f"{MODEL_DIR}/{MODEL_STEM}.{iso}-*.nc"):
        print(f"SKIP {ymd}: no model H1 file", flush=True); return

    # cam h1 files stamped at end of ouytput. FIlename carries the first sample time. 
    # the model steps that line up with a TEMPO/TROPOMI overpass get regridded   
    days = [(d0 + timedelta(days=k)).strftime("%Y-%m-%d") for k in (-1, 0, 1)]
    mfiles = sorted(f for day in days
                    for f in glob.glob(f"{MODEL_DIR}/{MODEL_STEM}.{day}-*.nc")) 
    cd = copy.deepcopy(yaml.safe_load(open(BASE)))
    cd["analysis"]["start_time"] = iso
    cd["analysis"]["end_time"] = f"{iso} 23:59:00"
    for k in ("output_dir", "output_dir_save", "output_dir_read"):
        cd["analysis"][k] = OUTDIR
    cd["analysis"]["save"]["paired"]["prefix"] = ymd
    cd["model"]["cam-chem-se"]["files"] = mfiles
    cd["model"]["cam-chem-se"]["scrip_file"] = _R["scrip"]

    # break job up 
    _grp = os.environ.get("OBS_GROUP", "").strip()
    keep = {s.strip() for s in _grp.split(",") if s.strip()} if _grp else None

    # vary regrid target 
    _rt = os.environ.get("REGRID_TARGET", "").strip()
    rt_val = ([s.strip() for s in _rt.split(",")] if "," in _rt else _rt) if _rt else None

    # vary regrid method (e.g. REGRID_METHOD=conservative to test the conservative
    # regridder over the full CONUS grid without editing the control)
    _rm = os.environ.get("REGRID_METHOD", "").strip()
    
    mapping = cd["model"]["cam-chem-se"].get("mapping", {})
    present = []
    for obs_name, src in OBS_SOURCES.items():
        if keep is not None and obs_name not in keep:
            cd["obs"].pop(obs_name, None); mapping.pop(obs_name, None)
            continue
        pat = f"{src['dir']}/{src['glob'].format(ymd=ymd)}"
        if obs_name in cd["obs"] and glob.glob(pat):
            cd["obs"][obs_name]["filename"] = pat
            if rt_val is not None:
                cd["obs"][obs_name]["regrid_target"] = rt_val
            if _rm:
                cd["obs"][obs_name]["regrid_method"] = _rm
                            
            present.append(obs_name)
        else:
            cd["obs"].pop(obs_name, None); mapping.pop(obs_name, None)
            print(f"  {ymd}: no {obs_name} files -> dropped", flush=True)
    if not present:
        print(f"SKIP {ymd}: no obs products with data", flush=True); return

    tmp = f"{YAMLDIR}/control_hcho_{ymd}.yaml"
    yaml.safe_dump(cd, open(tmp, "w"), sort_keys=False)
    print(f"==== pairing {iso}: {len(mfiles)} model | {present} ====", flush=True)
    
    an = driver.analysis()
    an.control = tmp
    an.read_control()
    an.open_models()
    an.open_obs()
    an.pair_data()
    an.save_analysis()
    
    print(f"==== DONE {iso} in {time.time()-t0:.0f}s ====", flush=True)
    for label, p in an.paired.items():
        print(f"  {label}: {dict(p.obj.sizes)}", flush=True)

if __name__ == "__main__":
    main()
