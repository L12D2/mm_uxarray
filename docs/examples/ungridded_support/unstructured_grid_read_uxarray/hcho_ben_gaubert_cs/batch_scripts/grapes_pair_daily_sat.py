"""
pair ONE day of TEMPO L2 HCHO vs CESM-SE. Driven by env YMD
"""

import os, glob, copy, time
from datetime import datetime
import yaml
from melodies_monet import driver

BASE = ("/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/satellite_grapes/control_grapes.yaml")

OUTDIR = ("/glade/work/lcthompson/mm_output/grapes")

YAMLDIR = OUTDIR + "/hcho_yaml"

MODEL_DIR = ("/glade/campaign/acom/acom-da/conus_outputs/"
             "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.nox_grapes.001/H1")

MODEL_STEM = "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.nox_grapes.001.cam.h1"

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

}

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(YAMLDIR, exist_ok=True)

def main():
    ymd = os.environ["YMD"]
    iso = datetime.strptime(ymd, "%Y%m%d").strftime("%Y-%m-%d")
    t0 = time.time()
    mfiles = sorted(glob.glob(f"{MODEL_DIR}/{MODEL_STEM}.{iso}-*.nc"))
    if not mfiles:
        print(f"SKIP {ymd}: no model H1 file", flush=True); return

    cd = copy.deepcopy(yaml.safe_load(open(BASE)))
    cd["analysis"]["start_time"] = iso
    cd["analysis"]["end_time"] = f"{iso} 23:59:00"
    for k in ("output_dir", "output_dir_save", "output_dir_read"):
        cd["analysis"][k] = OUTDIR
    cd["analysis"]["save"]["paired"]["prefix"] = ymd
    cd["model"]["cam-chem-se"]["files"] = mfiles

    # break job up 
    _grp = os.environ.get("OBS_GROUP", "").strip()
    keep = {s.strip() for s in _grp.split(",") if s.strip()} if _grp else None

    # vary regrid target 
    _rt = os.environ.get("REGRID_TARGET", "").strip()
    rt_val = ([s.strip() for s in _rt.split(",")] if "," in _rt else _rt) if _rt else None
    
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
