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

# CONTROL = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/control_tempo_l2_hcho_cesm_se.yaml"

CONTROL = os.environ.get("MM_CONTROL",
    "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/batch_scripts/control_master.yaml",)

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

def _swath_vars(p):
    """(model_var, obs_var, unc_var) for a native-swath pair."""
    dvs = list(p.obj.data_vars)
    unc = next((v for v in dvs if v.endswith(("_uncertainty", "_precision"))), None)
    mvs = getattr(p, "model_vars", None) or []
    model = next((v for v in mvs if v in dvs), None)
    if model is None:
        model = next((v for v in dvs if v != unc and "column" not in v), None)
    obs = next((v for v in dvs if v not in (model, unc)), None)
    return model, obs, unc


def _plot_swath(an, t0):
    """Standalone native-TEMPO swath plots (scatter + oversampling density).

    Triggered by env PLOT_SWATH=1. Iterates the *_swath pairs in an.paired and
    renders, per pair, a 3-panel obs/model/bias pixel scatter and a pixel-count
    oversampling density map. Bypasses an.plotting() -- the swath pixel vector
    (dims 'obs') does not fit the x/time spatial dispatch.
    """
    from melodies_monet.plots import satplots as splots

    out = an.control_dict.get("analysis", {}).get("output_dir") or "."
    os.makedirs(out, exist_ok=True)
    n = 0
    for label, p in an.paired.items():
        if not label.endswith("_swath") or "obs" not in getattr(p.obj, "dims", {}):
            continue
        model_var, obs_var, unc_var = _swath_vars(p)
        if model_var is None or obs_var is None:
            print(f"[swath] {label}: could not resolve model/obs vars "
                  f"from {list(p.obj.data_vars)}; skipping", flush=True)
            continue
        ylabel = p.obj[obs_var].attrs.get("units", "")
        print(f"[swath] {label}: model={model_var} obs={obs_var} "
              f"unc={unc_var} n={p.obj.sizes.get('obs')}", flush=True)
        splots.plot_swath_scatter(
            p.obj, model_var, obs_var, unc_var=unc_var,
            label_m=p.model, label_o=p.obs, ylabel=ylabel,
            outname=f"{out}/{label}.scatter",
        )
        splots.plot_swath_oversampling(
            p.obj, outname=f"{out}/{label}.oversampling",
        )
        n += 1
    print(f"[swath] plotted {n} swath pair(s) in {time.time()-t0:.0f}s", flush=True)


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

    _od = cd.get("analysis", {}).get("output_dir")
    if _od:
        os.makedirs(_od, exist_ok=True)
        
    print(f"[plot] opening model... ({time.time()-t0:.0f}s)", flush=True)
    an.open_models()
    print(f"[plot] reading paired... ({time.time()-t0:.0f}s)", flush=True)
    an.read_analysis()

    # inject different nicknames so we can get around the p.model 
    _nick = (an.control_dict.get("analysis", {}) or {}).get("sim_nickname")
    if _nick and len(an.models) == 1:
        _k = next(iter(an.models))
        an.models = {_nick: an.models[_k]}
        an.models[_nick].label = _nick
        _md = an.control_dict.get("model")
        if isinstance(_md, dict) and len(_md) == 1:
            an.control_dict["model"] = {_nick: next(iter(_md.values()))}
        for _p in an.paired.values():
            _p.model = _nick
        print(f"[plot] model label -> {_nick!r}", flush=True)

    if not an.paired:
        # Bail cleanly instead of an.plotting() indexing an empty pair list.
        print("[plot] no paired data after grid/only filter; nothing to plot.",
              flush=True)
        return
        
    if os.environ.get("PLOT_SWATH", "").strip():
        _plot_swath(an, t0)
        print(f"[plot] DONE (swath) in {time.time()-t0:.0f}s", flush=True)
        return

    for lbl, p in an.paired.items():
        print("Clipping...")
        def clip(var, lo, hi):
            if var not in p.obj:
                return
            da = p.obj[var]
            a = np.asarray(da.values)                      # materialize once
            np.putmask(a, (a <= lo) | (a >= hi), np.nan)    # mask in place ; no float duplicate
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
    #     qsub -N p_${tag} -o p_${tag}.log -l select=1:ncpus=1:mem=48GB -V submit_plot_sat.sh
    #   done
    # done

    

    # # gen controls
    # python make_native_controls.py                           
    
    # # plot
    # for RUN in nonbiog biog grapes mxcat; do
    #   export MM_CONTROL=$PWD/control_${RUN}_native.yaml PLOT_GRID=obs PLOT_ONLY=grp_native
    #   qsub -N pnat_${RUN} -o pnat_${RUN}.log -l select=1:ncpus=1:mem=64GB -V submit_plot_sat.sh
    # done

    
    # python make_controls.py                     
    # for RUN in nonbiog biog grapes mxcat; do
    #   for ONLY in grp5 grp6 grp7; do
    #     export MM_CONTROL=$PWD/control_${RUN}.yaml PLOT_GRID=model PLOT_ONLY=$ONLY
    #     tag=${RUN}_model_${ONLY}
    #     qsub -N p_$tag -o p_$tag.log -l select=1:ncpus=1:mem=64GB -V submit_plot_sat.sh
    #   done
    # done


    # # ### addtl regridding options 
    
    # # mxcat model-space
    # qsub -N pm_mxcat_cons -o pm_mxcat_cons.log \
    #      -l select=1:ncpus=1:mem=200GB -l walltime=24:00:00 \
    #      -v RUN=mxcat,REGRID_METHOD=conservative,REGRID_TARGET=model \
    #      submit_pair_tempo_native.sh
    
    # # ne0CONUS model-space
    # for RUN in nonbiog biog grapes; do
    #   qsub -N pm_${RUN}_cons -o pm_${RUN}_cons.log \
    #        -l select=1:ncpus=1:mem=128GB -l walltime=24:00:00 \
    #        -v RUN=$RUN,REGRID_METHOD=conservative,REGRID_TARGET=model \
    #        submit_pair_tempo_native.sh
    # done
    
    # # obs-grid native, conservative — CONUS cities on the ne0CONUS runs
    # for RUN in nonbiog biog grapes; do
    #   for CITY in atl dfw la den; do
    #     qsub -N po_${RUN}_${CITY}_cons -o po_${RUN}_${CITY}_cons.log \
    #          -l select=1:ncpus=1:mem=96GB -l walltime=24:00:00 \
    #          -v RUN=$RUN,CITY=$CITY,REGRID_METHOD=conservative,REGRID_TARGET=obs \
    #          submit_pair_tempo_native.sh
    #   done
    # done
    # # obs-grid native, conservative — mxcat over ATL + MEX
    # for CITY in atl mex; do
    #   qsub -N po_mxcat_${CITY}_cons -o po_mxcat_${CITY}_cons.log \
    #        -l select=1:ncpus=1:mem=128GB -l walltime=24:00:00 \
    #        -v RUN=mxcat,CITY=$CITY,REGRID_METHOD=conservative,REGRID_TARGET=obs \
    #        submit_pair_tempo_native.sh
    # done

