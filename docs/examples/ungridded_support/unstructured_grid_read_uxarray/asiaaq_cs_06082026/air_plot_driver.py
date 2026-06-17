import os, re, glob, copy, yaml, json, shutil, xarray as xr
import netCDF4
from datetime import datetime
from melodies_monet import driver

import sys
sys.path.insert(0, "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/output")
from montage_dc8 import make_montages

CTRL    = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/control_air_full_camp_cesm-se.yaml"
PAIRDIR = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/output/dc8_full_camp_pair"
PLOTROOT = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/output/air_full_plot"
TMPDIR = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/output/dc8_plot_yaml"

os.makedirs(TMPDIR, exist_ok=True)

base = yaml.safe_load(open(CTRL))

# flights that produced paired files 
dates = sorted(re.search(r'asiaaq_(\d{8})_dc8_', os.path.basename(f)).group(1)
               for f in glob.glob(f"{PAIRDIR}/asiaaq_*_dc8_cam-chem-se-era5.nc4"))
print("flights to plot:", dates)

for ymd in dates:
    iso = datetime.strptime(ymd, "%Y%m%d").strftime("%Y-%m-%d")
    cd  = copy.deepcopy(base)

    outdir = f"{PLOTROOT}/{ymd}"
    os.makedirs(outdir, exist_ok=True)
    cd["analysis"]["output_dir"]      = outdir
    cd["analysis"]["output_dir_save"] = outdir
    cd["analysis"]["output_dir_read"] = PAIRDIR
    cd["analysis"]["start_time"]      = iso
    cd["analysis"]["end_time"]        = f"{iso} 23:59:00"

    # this flight's two paired files only
    fn = cd["analysis"]["read"]["paired"]["filenames"]
    for model in ["cam-chem-se-era5", "cam-chem-se-merra2"]:
        f = f"{PAIRDIR}/asiaaq_{ymd}_dc8_{model}.nc4"
        fn[f"dc8_{model}"] = [f] if os.path.exists(f) else []
    fn = {k: v for k, v in fn.items() if v}              # drop any model w/o a file
    cd["analysis"]["read"]["paired"]["filenames"] = fn
    avail = list(fn.keys())                              # pairs present for THIS flight

    # way to skip variables that were not recorded on flight
    MODELCOL = {"O3": "O3_new", "CO": "CO_new", "NO2": "NO2_new", "temperature": "T"}   # obs_var to model column
    COORDS   = {"pressure_obs", "altitude"}                         # always keep (vertical axis)

    present = set(MODELCOL)
    for files in fn.values():                       # every available pair file for this flight
        ds = xr.open_dataset(files[0])
        dv = set(ds.data_vars)
        for ov, mcol in MODELCOL.items():
            if mcol not in dv or int(ds[mcol].count()) == 0:   # model col absent/all-NaN
                present.discard(ov)
        ds.close()
        
    print(f"   {ymd}: plotting vars -> {sorted(present)}", flush=True)
    if not present:
        print(f"   SKIP {ymd}: nothing to plot"); continue

    PRUNEDIR = f"{TMPDIR}/pruned"
    os.makedirs(PRUNEDIR, exist_ok=True)
    new_fn = {}
    
    for key, files in fn.items():
        dst = f"{PRUNEDIR}/{ymd}_{os.path.basename(files[0])}"
        shutil.copyfile(files[0], dst)
        with netCDF4.Dataset(dst, "a") as nc:
            meta = json.loads(nc.getncattr("dict_json"))
            ov, mv = meta["obs_vars"], meta["model_vars"]
            idx = [i for i, o in enumerate(ov) if o in present]   # keep aligned pairs
            meta["obs_vars"]   = [ov[i] for i in idx]
            meta["model_vars"] = [mv[i] for i in idx]
            nc.setncattr("dict_json", json.dumps(meta))
        new_fn[key] = [dst]
    cd["analysis"]["read"]["paired"]["filenames"] = new_fn
    
    vd = cd["obs"]["dc8"]["variables"]
    cd["obs"]["dc8"]["variables"] = {k: v for k, v in vd.items() if k in COORDS or k in present}
    print(f"   {ymd}: plotting vars -> {sorted(present)}", flush=True)
    if not present:
        print(f"   SKIP {ymd}: no gas variables present"); continue
        
    # label every plot/stat with the date, and only reference available pairs
    for g in cd["plots"].values():
        g["domain_name"] = [f"DC8 {iso}"]
        g["data"] = [d for d in g["data"] if d in avail]
    cd["stats"]["domain_name"] = [f"DC8 {iso}"]
    cd["stats"]["data"] = [d for d in cd["stats"]["data"] if d in avail]

    tmp = f"{TMPDIR}/control_dc8_plot_{ymd}.yaml"
    yaml.safe_dump(cd, open(tmp, "w"), sort_keys=False)

    print(f"==== plotting flight {iso} ====", flush=True)
    an = driver.analysis()
    an.control = tmp
    an.read_control()
    an.read_analysis()
    an.plotting()
    an.stats()

make_montages(PLOTROOT, f"{os.path.dirname(PLOTROOT)}/dc8_montages")