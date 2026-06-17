import glob, os, re, copy, yaml, xarray as xr, pandas as pd
from datetime import datetime, timedelta
from melodies_monet import driver

BASE = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/write_met_data/full_camp/dc8_pair.yaml"   

# make directory for the yamls 
os.makedirs("/glade/derecho/scratch/lcthompson", exist_ok=True)

# make output dir from the yaml 
os.makedirs("/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/output/dc8_full_camp_pair", exist_ok=True)


ICT  = "/glade/campaign/acom/acom-weather/emmons/ASIAAQ_obs/DC8"

ERA5_HIST = "/glade/campaign/acom/acom-weather/emmons/ASIAAQ_sims/f.e3b08.FHISTC_LTt1s.camsm11_finn27.era5.GEMSne30x8.2024.01/atm/hist"
ERA5_STEM = "f.e3b08.FHISTC_LTt1s.camsm11_finn27.era5.GEMSne30x8.2024.01.cam.h5a"

MER_HIST  = "/glade/campaign/acom/acom-weather/emmons/ASIAAQ_sims/f.e3b08.FHISTC_LTt1s.camsm11_finn27.2024.merra2.ne30.01/atm/hist"
MER_STEM  = "f.e3b08.FHISTC_LTt1s.camsm11_finn27.2024.merra2.ne30.01.cam.h4i"
#OBSDIR    = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/preprocessing/dc8_data"

TMPDIR = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/output/dc8_yaml"
os.makedirs(TMPDIR, exist_ok=True)

def mfile(hist, stem, ymd):
    return f"{hist}/{stem}.{ymd[:4]}-{ymd[4:6]}-{ymd[6:8]}-03600.nc"

# flight dates straight from the .ict filenames
dates = sorted({re.search(r'_(\d{8})_', os.path.basename(f)).group(1)
                for f in glob.glob(f"{ICT}/asiaaq-mrg10_dc8_*.ict")})

# debug small patch 
#dates = ['20240213']

print("flights:", dates)

base = yaml.safe_load(open(BASE))

for ymd in dates:
    d   = datetime.strptime(ymd, "%Y%m%d")
    nxt = (d + timedelta(days=1)).strftime("%Y%m%d")
    iso = d.strftime("%Y-%m-%d")

    cd = copy.deepcopy(base)
    cd["analysis"]["start_time"] = iso
    cd["analysis"]["end_time"]   = f"{iso} 23:59:00"
    cd["analysis"]["save"]["paired"]["prefix"] = f"asiaaq_{ymd}"

    cd["model"]["cam-chem-se-era5"]["files"]   = [mfile(ERA5_HIST, ERA5_STEM, x) for x in (ymd, nxt)]
    cd["model"]["cam-chem-se-merra2"]["files"] = [mfile(MER_HIST,  MER_STEM,  x) for x in (ymd, nxt)]

    # point obs at the .ict since mm has an ict reader
    ict = glob.glob(f"{ICT}/asiaaq-mrg10_dc8_{ymd}_*.ict")
    if not ict:
        print(f"  SKIP {ymd}: no .ict"); continue
    cd["obs"]["dc8"]["filename"] = ict[0]

    for m in cd["model"].values():
        m["files"] = [p for p in m["files"] if os.path.exists(p)]
    if not all(m["files"] for m in cd["model"].values()):
        print(f"  SKIP {ymd}: missing model file"); continue

    tmp = f"{TMPDIR}/control_dc8_{ymd}.yaml"
    yaml.safe_dump(cd, open(tmp, "w"), sort_keys=False)

    print(f"==== pairing flight {iso} ====")
    an = driver.analysis()
    an.control = tmp
    an.read_control()
    an.open_models()
    an.open_obs()
    an.pair_data()
    an.save_analysis() 

    outdir  = cd["analysis"]["output_dir_save"]
    
    written = sorted(glob.glob(f"{outdir}/asiaaq_{ymd}_*"))
    print(f"==== DONE {iso}: {len(written)} paired file(s) -> {outdir}", flush=True)
    for w in written:
        print("       ", os.path.basename(w), flush=True)
