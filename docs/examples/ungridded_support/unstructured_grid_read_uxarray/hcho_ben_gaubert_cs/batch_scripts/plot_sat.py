"""

python pair_sat.py

Produces in output_dir_save:
    
"""
import time
import warnings
import numpy as np

warnings.filterwarnings("ignore")

from melodies_monet import driver

CONTROL = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/control_tempo_l2_hcho_cesm_se.yaml"

def main():
    t0 = time.time()
    an = driver.analysis()
    an.control = CONTROL
    an.read_control()

    print("[pair] opening models (lazy)...", flush=True)
    an.open_models()
    
    # print(f"[pair] opening obs... ({time.time() - t0:.0f}s)", flush=True)
    # an.open_obs()

    print(f"[reading] reading paired nc4... ({time.time() - t0:.0f}s)", flush=True)
    #an.pair_data()
    an.read_analysis()

    # try and reduce some noise 
    for lbl, p in an.paired.items():
        for v in ("vertical_column", "formaldehyde_tropospheric_vertical_column", "CH2O"):
            if v in p.obj:
                vc = p.obj[v]
                p.obj[v] = vc.where(np.isfinite(vc) & (vc > -2e16) & (vc < 5e16))
    
    print(f"[saved pair] plotting... ({time.time() - t0:.0f}s)", flush=True)
    an.plotting()

    an.stats()

    print(f"[plotting] DONE in {time.time() - t0:.0f}s", flush=True)
    for label, p in an.paired.items():
        print(f"  {label}: {dict(p.obj.sizes)}", flush=True)

if __name__ == "__main__":
    main()
