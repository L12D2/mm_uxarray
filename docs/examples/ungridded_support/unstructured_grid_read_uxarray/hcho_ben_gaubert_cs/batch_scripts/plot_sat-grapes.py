"""

python pair_sat.py

Produces in output_dir_save:
    
"""
import time
import os
import warnings
import numpy as np

warnings.filterwarnings("ignore")

from melodies_monet import driver

CONTROL = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/satellite_grapes/control_grapes.yaml"

def _grid_of(label):
    """Which grid a paired label belongs to."""
    if label.endswith("_obsgrid") or label.endswith("_series"):
        return "obs"        # lat/lon L3 + weighted series 
    return "model"          # unstruct mesh

def _keep_plot(name, g, grid, only):
    data = g.get("data") or []
    if not data:
        return False
    if grid in ("obs", "model") and not all(_grid_of(d) == grid for d in data):
        return False
    if only and only not in name:
        return False
    return True

def main():
    t0 = time.time()
    grid = os.environ.get("PLOT_GRID", "").strip()   # "obs" | "model" | ""
    only = os.environ.get("PLOT_ONLY", "").strip()   # group-name substring, e.g. "grp4_"

    an = driver.analysis()
    an.control = CONTROL
    an.read_control()

    cd = an.control_dict
    
    if grid or only:
        cd["plots"] = {n: g for n, g in cd.get("plots", {}).items()
                       if _keep_plot(n, g, grid, only)}
        # read ONLY the labels the surviving plots need -> minimal memory
        used = set()
        for g in cd["plots"].values():
            used.update(g["data"])
        fn = cd["analysis"]["read"]["paired"]["filenames"]
        cd["analysis"]["read"]["paired"]["filenames"] = {k: v for k, v in fn.items() if k in used}
        if cd.get("stats", {}).get("data"):
            cd["stats"]["data"] = [d for d in cd["stats"]["data"] if d in used]
        print(f"[plot] grid={grid!r} only={only!r}: {len(used)} labels, "
              f"{len(cd['plots'])} plot groups", flush=True)

    print(f"[plot] opening model... ({time.time()-t0:.0f}s)", flush=True)
    an.open_models()
    print(f"[plot] reading paired... ({time.time()-t0:.0f}s)", flush=True)
    an.read_analysis()

    for lbl, p in an.paired.items():
        print("Clipping...")
        def clip(var, lo, hi):
            if var not in p.obj:
                return
            da = p.obj[var]
            a = np.asarray(da.values)                      # materialize once
            np.putmask(a, (a <= lo) | (a >= hi), np.nan)    # mask in place -- no float duplicate
            p.obj[var] = da.copy(data=a)                    
        for v in ("vertical_column", "formaldehyde_tropospheric_vertical_column", "CH2O"):
            clip(v, -2e16, 5e16)
        for v in ("vertical_column_troposphere", "nitrogendioxide_tropospheric_column", "NO2"):
            clip(v, -2e15, 5e16)
        for v in ("carbonmonoxide_total_column", "CO"):
            clip(v, 0, 1e19)
            
    print(f"[plot] plotting... ({time.time()-t0:.0f}s)", flush=True)
    an.plotting()
    print(f"[plot] DONE in {time.time()-t0:.0f}s", flush=True)
    for label, p in an.paired.items():
        print(f"  {label}: {dict(p.obj.sizes)}", flush=True)


if __name__ == "__main__":
    main()


    # for lbl, p in an.paired.items():
    #     for v in ("vertical_column", "formaldehyde_tropospheric_vertical_column", "CH2O"):
    #         clip_and_report(p, lbl, v, -2e16, 5e16)   # HCHO
    #     for v in ("vertical_column_troposphere", "nitrogendioxide_tropospheric_column", "NO2"):
    #         clip_and_report(p, lbl, v, -2e15, 5e16)   # NO2
    #     for v in ("carbonmonoxide_total_column", "CO"):
    #         clip_and_report(p, lbl, v, 0, 1e19)       # CO
    
    # for GRID in obs model; do
    #   for ONLY in grp1_ grp1b_ grp2_ grp3 grp4_ grp5; do
    #     export PLOT_GRID=$GRID PLOT_ONLY=$ONLY
    #     tag=${GRID}_${ONLY%_}
    #     qsub -N p_${tag} -o p_${tag}.log -l select=1:ncpus=1:mem=48GB -V submit_plot_sat-grapes.sh
    #   done
    # done

# qsub -N p_co -o p_co.log -l select=1:ncpus=1:mem=48GB \
#      -v PLOT_ONLY=_co submit_plot_sat.sh
