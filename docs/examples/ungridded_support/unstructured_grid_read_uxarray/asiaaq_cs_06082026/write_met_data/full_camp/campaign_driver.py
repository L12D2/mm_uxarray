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
OBSDIR    = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/preprocessing/dc8_data"

# need to split the merge all .nc 
# this file DOES NOT read .ict files 

ds = xr.open_dataset(f"{OBSDIR}/asiaaq_dc8_merge_all.nc")
t  = pd.to_datetime(ds["time"].values)
for day, idx in pd.Series(range(len(t)), index=t).groupby(t.normalize()):
    ds.isel(time=idx.values).to_netcdf(f"{OBSDIR}/asiaaq_dc8_{day:%Y%m%d}.nc")
    print("wrote", f"asiaaq_dc8_{day:%Y%m%d}.nc", len(idx))


def mfile(hist, stem, ymd):
    return f"{hist}/{stem}.{ymd[:4]}-{ymd[4:6]}-{ymd[6:8]}-03600.nc"

# flight dates straight from the .ict filenames
dates = sorted({re.search(r'_(\d{8})_', os.path.basename(f)).group(1)
                for f in glob.glob(f"{ICT}/asiaaq-mrg10_dc8_*.ict")})

# debug small patch 
# dates = [d for d in sorted({re.search(r'_(\d{8})_', os.path.basename(f)).group(1)
#                             for f in glob.glob(f"{ICT}/asiaaq-mrg10_dc8_*.ict")})
#          if d < "20240310"]

print("flights:", dates)

base = yaml.safe_load(open(BASE))

for ymd in dates:
    d   = datetime.strptime(ymd, "%Y%m%d")
    nxt = (d + timedelta(days=1)).strftime("%Y%m%d")   # cover midnight-crossing flights
    iso = d.strftime("%Y-%m-%d")

    cd = copy.deepcopy(base)
    cd["analysis"]["start_time"] = iso
    cd["analysis"]["end_time"]   = f"{iso} 23:59:00"
    cd["analysis"]["save"]["paired"]["prefix"] = f"asiaaq_{ymd}"

    cd["model"]["cam-chem-se-era5"]["files"]   = [mfile(ERA5_HIST, ERA5_STEM, x) for x in (ymd, nxt)]
    cd["model"]["cam-chem-se-merra2"]["files"] = [mfile(MER_HIST,  MER_STEM,  x) for x in (ymd, nxt)]

    # per-flight obs file
    perflight = f"{OBSDIR}/asiaaq_dc8_{ymd}.nc"
    cd["obs"]["dc8"]["filename"] = perflight if os.path.exists(perflight) \
        else f"{OBSDIR}/asiaaq_dc8_merge_all.nc"

    # keep only model files that actually exist
    for m in cd["model"].values():
        m["files"] = [p for p in m["files"] if os.path.exists(p)]
    if not all(m["files"] for m in cd["model"].values()):
        print(f"  SKIP {ymd}: missing model file"); continue

    tmp = f"/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/output/dc8_yaml/control_dc8_{ymd}.yaml"
    yaml.safe_dump(cd, open(tmp, "w"), sort_keys=False)

    print(f"==== pairing flight {iso} ====")
    an = driver.analysis()
    an.control = tmp
    an.read_control()
    an.open_models()
    an.open_obs()
    an.pair_data()
    an.save_analysis()      
