"""

python pair_aircraft.py

Produces in output_dir_save:
    asiaaq_dc8_dc8_cam-chem-se-era5.nc4
    asiaaq_dc8_dc8_cam-chem-se-merra2.nc4
    
"""
import time
import warnings

warnings.filterwarnings("ignore")

from melodies_monet import driver

CONTROL = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/write_met_data/batch_scripts/sfc_air_pair.yaml"

def main():
    t0 = time.time()
    an = driver.analysis()
    an.control = CONTROL
    an.read_control()

    print("[pair] opening models (lazy)...", flush=True)
    an.open_models()
    print(f"[pair] opening obs... ({time.time() - t0:.0f}s)", flush=True)
        
    an.open_obs()

    print(f"[pair] pairing... ({time.time() - t0:.0f}s)", flush=True)
    an.pair_data()

    print(f"[pair] saving paired NetCDFs... ({time.time() - t0:.0f}s)", flush=True)
    an.save_analysis()

    print(f"[pair] DONE in {time.time() - t0:.0f}s", flush=True)
    for label, p in an.paired.items():
        print(f"  {label}: {dict(p.obj.sizes)}", flush=True)

if __name__ == "__main__":
    main()
