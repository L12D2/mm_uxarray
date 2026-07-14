from melodies_monet import driver
import os, re, glob
from datetime import datetime, timedelta


def extend_model_files_prev_day(control_dict, ndays=1):
    """Prepend the model file(s) for the day(s) BEFORE analysis.start_time.

    CAM h1 files are stamped at the END of the output interval, and the
    filename carries the FIRST sample time -- so a given day's 00 UTC record
    lives in the PREVIOUS day's file (its last record) 
    """
    mblock = next(iter(control_dict["model"].values()))
    fstr = mblock.get("files")
    if not isinstance(fstr, str):
        return                                      # already explicit list; leave as-is
    files = set(glob.glob(os.path.expandvars(fstr)))
    start = str(control_dict["analysis"]["start_time"])[:10]     # 'YYYY-MM-DD'
    d0 = datetime.strptime(start, "%Y-%m-%d")
    d, b = os.path.dirname(fstr), os.path.basename(fstr)
    m = re.match(r"(.+\.cam\.h\d(?:i)?)\.", b)       # capture 'STEM.cam.h1' / '...h3i'
    if m:
        prefix = m.group(1)
        for k in range(1, ndays + 1):
            prev = (d0 - timedelta(days=k)).strftime("%Y-%m-%d")
            files |= set(glob.glob(os.path.join(d, f"{prefix}.{prev}-*.nc")))
    mblock["files"] = sorted(files)                  # list -> used as literal paths
    
BASE = ("/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/"
        "unstructured_grid_read_uxarray/hcho_ben_gaubert_cs")

CONTROLS = {
    "nonbiog": f"{BASE}/sfc_img/control_airnow_cesm_se.yaml",
    "biog":    f"{BASE}/sfc_bio/control_airnow_cesm_se-bio.yaml",
    "grapes":  f"{BASE}/sfc_grapes/grapes_control_airnow_cesm_se.yaml",
    "mxcat":   f"{BASE}/sfc_mxcat/mxcat_control_airnow_cesm_se.yaml",
}

RUN = os.environ.get("RUN", "nonbiog").strip()
if RUN not in CONTROLS:
    raise SystemExit(f"RUN={RUN!r} not in {sorted(CONTROLS)}")

def main():
    an = driver.analysis()
    an.control = CONTROLS[RUN]
    print(f"==== sfc pairing RUN={RUN}: {an.control} ====", flush=True)

    an.read_control()
    
    # 00 UTC fix: pull in the previous day's model file so the first day's 00Z
    # has a model value to pair against (CAM end-of-interval stamping).
    extend_model_files_prev_day(an.control_dict)
    
    an.open_models()
    an.open_obs()

    print("pairing data....", flush=True)
    an.pair_data()
    an.save_analysis()

    # print("plotting data....", flush=True)
    # an.plotting()
    print("done!")


if __name__ == "__main__":
    main()
    

# def main():
#     # run mm
#     an = driver.analysis()
#     an.control = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/sfc_img/control_airnow_cesm_se.yaml"

#     an.control = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/sfc_mxcat/mxcat_control_airnow_cesm_se.yaml"

#     an.control = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/sfc_bio/control_airnow_cesm_se-bio.yaml"

#     an.control = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/sfc_grapes/grapes_control_airnow_cesm_se.yaml"
    
#     an.read_control()
#     an.open_models()
#     an.open_obs()
    
#     print("pairing data....", flush=True)
#     an.pair_data()
#     an.save_analysis()
#     #an.read_analysis()
    
#     # print("plotting data....", flush=True)
#     an.plotting()
#     # an.stats()

# if __name__ == "__main__":
#     main()
