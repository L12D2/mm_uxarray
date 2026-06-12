
"""

Batch run: TEMPO L2 NO2 vs CESM-fv (CONUS)

"""
import time
import warnings

warnings.filterwarnings("ignore")

from melodies_monet import driver

CONTROL = "./control_tempo_l2_no2_cesm_se.yaml"


def main():
    t0 = time.time()
    an = driver.analysis()
    an.control = CONTROL
    an.read_control()

    print("[tempo] open_models (lazy)...", flush=True)
    an.open_models()
    for m in an.models.values():
        if "lon" in m.obj.variables and "longitude" not in m.obj.variables:
            m.obj = m.obj.rename({"lon": "longitude", "lat": "latitude"})
        
    print(f"[tempo] open_obs... ({time.time() - t0:.0f}s)", flush=True)
    an.open_obs()

    print(f"[tempo] pair_data (swath regrid -> unstructured)... ({time.time() - t0:.0f}s)", flush=True)
    an.pair_data()

    print(f"[tempo] save_analysis... ({time.time() - t0:.0f}s)", flush=True)
    an.save_analysis()                      # paired saved -- safe past this point

    for label, p in an.paired.items():
        print(f"  paired {label}: {dict(p.obj.sizes)}", flush=True)

    try:
        print(f"[tempo] plotting... ({time.time() - t0:.0f}s)", flush=True)
        an.plotting()
    except Exception as e:
        print(f"[tempo] plotting FAILED ({type(e).__name__}): {e}", flush=True)

    try:
        print(f"[tempo] stats... ({time.time() - t0:.0f}s)", flush=True)
        an.stats()
    except Exception as e:
        print(f"[tempo] stats FAILED ({type(e).__name__}): {e}", flush=True)

    print(f"[tempo] DONE in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
