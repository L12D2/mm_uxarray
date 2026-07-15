"""Render per-run MELODIES-MONET control YAMLs from one shared master yaml

Help unify multiple yaml files if you are running the same model and obs with different emissions scenarios. 
Makes yaml maintenance easier. 

Usage:
    python make_controls.py            # writes control_<run>.yaml per run
    python make_controls.py mxcat      # just one run

The master (control_master.yaml) holds everything shared: obs blocks,
mapping, all plot groups, stats. runs.yaml holds the handful of per-run
differences: output dirs, model files, SCRIP, extra zoom domains.
Edit those two; NEVER edit the generated control_<run>.yaml files.

Per-run edits applied to a deep copy of the master:
  - analysis.output_dir / output_dir_save / output_dir_read, montage dirs
  - read.paired.filenames: every label repointed at <paired_dir>/<prefix>*_<label>.nc4
  - model files + scrip_file (first model entry)
  - add_domains: appended to every plot group (and stats) that already uses
    an auto-region:box domain, so zoom cities show up everywhere a box does.
"""

# run 

# ./hcho_ben_gaubert_cs/batch_scripts
# python make_controls.py                      # control_{nonbiog,biog,grapes,mxcat}.yaml
# for RUN in nonbiog biog grapes mxcat; do
#   for ONLY in grp5 grp6 grp7; do
#     export MM_CONTROL=$PWD/control_${RUN}.yaml PLOT_GRID=model PLOT_ONLY=$ONLY
#     tag=${RUN}_model_${ONLY}
#     qsub -N p_$tag -o p_$tag.log -l select=1:ncpus=1:mem=48GB -V submit_plot_sat.sh
#   done
# done

# # cons mod space
# python make_controls.py nonbiog_cons biog_cons grapes_cons mxcat_cons
# for RUN in nonbiog_cons biog_cons grapes_cons mxcat_cons; do
#   for ONLY in grp6 grp7 grp8; do          # TEMPO model-grid groups ONLY
#     export MM_CONTROL=$PWD/control_${RUN}.yaml PLOT_GRID=model PLOT_ONLY=$ONLY
#     qsub -N pc_${RUN}_${ONLY} -o pc_${RUN}_${ONLY}.log -l select=1:ncpus=1:mem=40GB -V submit_plot_sat.sh
#   done
# done

# # native obs space
# python make_native_controls.py
# for RUN in nonbiog biog grapes mxcat; do
#   export MM_CONTROL=$PWD/control_${RUN}_native_cons.yaml PLOT_GRID=obs PLOT_ONLY=grp_native
#   qsub -N pnc_${RUN} -o pnc_${RUN}.log -l select=1:ncpus=1:mem=30GB -V submit_plot_sat.sh
# done
# # then flip METHOD="radius_mean" at top of make_native_controls.py, rerun it, and plot control_*_native_rmean.yaml

import copy
import pathlib
import sys
import os

import yaml

HERE = pathlib.Path(__file__).parent

METHOD_TAG = {"conservative": "cons", "conservative_normed": "consn",
              "radius_mean": "rmean", "nearest_s2d": "nn", "nearest_d2s": "nnd",
              "bilinear": "bilin"}
CITY_BOX = {"atl": "ATL_metro", "mex": "MEX_metro", "la": "LA_metro",
            "den": "DEN_metro", "dfw": "DFW_metro", "seus": "SE_US"}

class NoAliasDumper(yaml.SafeDumper):
    """Dump shared objects (from master YAML anchors) as plain copies, not *id aliases."""
    def ignore_aliases(self, data):
        return True

def render(master, name, spec, common=None):
    cd = copy.deepcopy(master)
    an = cd["analysis"]

    an["output_dir"] = spec["plot_dir"]
    an["output_dir_save"] = spec["paired_dir"]
    an["output_dir_read"] = spec["paired_dir"]
    if spec.get("nickname"):                      # display label for legends/titles
        an["sim_nickname"] = spec["nickname"]
    if "montage" in an:
        an["montage"]["plot_dir"] = spec["plot_dir"]
        an["montage"]["outdir"] = spec["plot_dir"] + "/montages"

    prefix = spec.get("paired_prefix") or an["save"]["paired"]["prefix"]
    fns = an["read"]["paired"]["filenames"]
    for label in fns:
        fns[label] = [f"{spec['paired_dir']}/{prefix}*_{label}.nc4"]

    mod = next(iter(cd["model"].values()))

    mod["files"] = spec["model_files"]
    mod["scrip_file"] = spec["scrip_file"]

    # instrument / target scoping 
    inst = spec.get("instrument")
    targets = spec.get("targets")

    def _label_ok(label):
        if inst and not label.startswith(f"{inst}_l2"):
            return False
        if targets is not None:
            tgt = ("obsgrid" if label.endswith("_obsgrid")
                   else "series" if label.endswith("_series") else "model")
            if tgt not in targets:
                return False
        return True

    if inst or targets is not None:
        an["read"]["paired"]["filenames"] = {
            lbl: v for lbl, v in fns.items() if _label_ok(lbl)}
        _plots = cd.get("plots", {})
        for gname in list(_plots):
            g = _plots[gname]
            g["data"] = [d for d in (g.get("data") or []) if _label_ok(d)]
            if not g["data"]:
                del _plots[gname]        # nothing left to plot for this run
        if isinstance(cd.get("stats"), dict):
            cd["stats"]["data"] = [d for d in (cd["stats"].get("data") or [])
                                   if _label_ok(d)]
            if not cd["stats"]["data"]:
                cd.pop("stats", None)
        print(f"  [{name}] scoped to instrument={inst!r} targets={targets!r}: "
              f"{len(an['read']['paired']['filenames'])} labels, "
              f"{len(cd.get('plots', {}))} plot groups")
        
    blocks = list(cd.get("plots", {}).values())

    if isinstance(cd.get("stats"), dict):
        blocks.append(cd["stats"])
    # shared domains first, then per-run (per-run keys win on conflict)
    add = dict((common or {}).get("add_domains") or {})
    add.update(spec.get("add_domains") or {})
    for dom, dspec in add.items():
        # dict with bounds -> box domain; null/other -> 'all' domain
        is_box = isinstance(dspec, dict) and "bounds" in dspec
        dtype = "auto-region:box" if is_box else "all"
        for g in blocks:
            types = g.get("domain_type") or []
            names = g.get("domain_name") or []
            # Augment a group if it already plots a spatial box, OR if its
            # domain_type is left empty (domains come entirely from add_domains).
            # Groups with a non-box domain set (e.g. ["all"]/CONUS series) are
            # left alone.
            if types and "auto-region:box" not in types:
                continue
            if dom in names:
                continue
            if is_box:
                if not isinstance(g.get("domain_info"), dict):
                    g["domain_info"] = {}
                g["domain_info"][dom] = dspec
            g["domain_type"] = list(types) + [dtype]
            g["domain_name"] = list(names) + [dom]

    # tidy: drop a leftover null domain_info; warn on any group that ended up
    # with no domains (e.g. emptied group but add_domains was empty)
    for gname, g in cd.get("plots", {}).items():
        if g.get("domain_info") is None:
            g.pop("domain_info", None)
        if not (g.get("domain_type") or []):
            print(f"  WARNING: plot group '{gname}' has no domains after "
                  "injection (empty domain_type + no add_domains).")

    out = HERE / f"control_{name}.yaml"
    header = (
        f"# GENERATED by make_controls.py (run: {name}) from "
        "control_master.yaml + runs.yaml.\n"
        "# Do not edit this file -- edit the master/runs and regenerate.\n"
    )
    with open(out, "w") as f:
        f.write(header)
        yaml.safe_dump(cd, f, sort_keys=False, width=120)
    print(f"wrote {out}")

def render_native(master, run, inst, product, ncfg, spec, boxes):
    """Render one native-resolution city control (obsgrid or swath).

    Same master skeleton as the standard controls, but the obs/model reads and
    plot groups are BUILT here (per city x species) for the native product,
    rather than reused from the master's authored plots. Emits
    control_{run}_[<inst>_]native_{mtag}[_swath].yaml.
    """
    mtag = METHOD_TAG.get(ncfg["method"], ncfg["method"][:5])
    icfg = ncfg["instruments"][inst]
    res = icfg["res"]
    restag = ("%g" % res).replace(".", "p")
    root = ncfg["root"]
    min_obs = ncfg.get("min_obs", 3)

    paired_dir = f"{root}/{run}_{inst}_native"
    plot_dir = spec["plot_dir"].rstrip("/") + (
        "_native" if inst == "tempo" else f"_{inst}_native")

    cd = copy.deepcopy(master)
    an = cd["analysis"]
    an["output_dir"] = plot_dir
    an["output_dir_save"] = paired_dir
    an["output_dir_read"] = paired_dir
    if spec.get("nickname"):                      # display label for legends/titles
        an["sim_nickname"] = spec["nickname"]
    if "montage" in an:
        an["montage"]["plot_dir"] = plot_dir
        an["montage"]["outdir"] = plot_dir + "/montages"
    an.pop("save", None)                       # plotting only

    mod = next(iter(cd["model"].values()))
    mod["files"] = spec["model_files"]
    mod["scrip_file"] = spec["scrip_file"]

    cd.pop("stats", None)                      # obs blocks kept (ylabel/vmin lookup)

    filenames, plots = {}, {}
    for city in ncfg["cities"].get(run, []):
        box = CITY_BOX.get(city)
        if box is None or box not in boxes:
            print(f"  WARNING: native {run}/{city}: box {box} not in add_domains; skipping.")
            continue
        for sp in icfg["species"]:
            obs_name = f"{inst}_l2_{sp}"
            label = f"{city}_{obs_name}_cam-chem-se_{product}"
            filenames[label] = [
                f"{paired_dir}/{city}_{restag}_{mtag}_*_{obs_name}_cam-chem-se_{product}.nc4"]
            common_grp = {
                "domain_type": ["auto-region:box"],
                "domain_name": [box],
                "domain_info": {box: boxes[box]},
                "data": [label],
            }
            if product == "swath":
                plots[f"grp_swath_{city}_{sp}"] = {
                    "type": "spatial_swath",
                    "fig_kwargs": {"figsize": [22, 6], "states": True},
                    "text_kwargs": {"fontsize": 14},
                    "data_proc": {"render": "auto"},
                    **common_grp,
                }
            else:                              # obsgrid -> gridded spatial_overlay
                plots[f"grp_native_{city}_{sp}"] = {
                    "type": "spatial_overlay",
                    "fig_kwargs": {"figsize": [18, 5], "cbar_orientation": "horizontal",
                                   "cbar_kwargs": {"shrink": 0.6},
                                   "states": True, "counties": True},
                    "text_kwargs": {"fontsize": 16},
                    "data_proc": {"time_reduction": "mean", "daily_first": True,
                                  "common_mask": True, "min_obs": min_obs,
                                  "set_axis": True, "rem_obs_nan": True},
                    **common_grp,
                }
    an["read"] = {"paired": {"method": "netcdf", "filenames": filenames}}
    cd["plots"] = plots

    itag = "" if inst == "tempo" else f"{inst}_"
    psfx = "" if product == "obsgrid" else f"_{product}"
    out = HERE / f"control_{run}_{itag}native_{mtag}{psfx}.yaml"
    with open(out, "w") as f:
        f.write(f"# GENERATED by make_controls.py (native {inst}/{product}, run: {run}).\n"
                "# Do not edit -- edit runs.yaml and regenerate.\n")
        yaml.safe_dump(cd, f, sort_keys=False, width=120)
    print(f"wrote {out}  ({len(plots)} groups, {len(filenames)} labels)")

def main():
    cfg = yaml.safe_load(open(HERE / "runs.yaml"))
    master = yaml.safe_load(open(HERE / cfg.get("master", "control_master.yaml")))
    common = cfg.get("common")
    boxes = (common or {}).get("add_domains") or {}
    args = sys.argv[1:]

    # standard (master-driven) controls
    wanted = args or list(cfg["runs"])

    nicknames = cfg.get("nicknames", {})
    def _nick(run):
        return nicknames.get(run.split("_")[0])
    for name in wanted:
        render(master, name, {**cfg["runs"][name], "nickname": _nick(name)},
               common=common)

    # surface controls
    sfc = cfg.get("sfc")
    if sfc and not args:
        sfc_master = yaml.safe_load(open(HERE / sfc.get("master", "control_sfc_master.yaml")))
        for run, sspec in sfc.get("runs", {}).items():
            base = cfg["runs"][run]
            spec = {
                "plot_dir": sspec["plot_dir"],
                "paired_dir": base["paired_dir"],
                "paired_prefix": sspec["prefix"],
                "model_files": base["model_files"],
                "scrip_file": base["scrip_file"],
                "nickname": _nick(run),
            }
            render(sfc_master, f"{run}_sfc", spec, common=common)
            
    # native city products -- only on a full run (no explicit run args)
    ncfg = cfg.get("native")
    if ncfg and not args:
        for run in ncfg.get("cities", {}):
            spec = {**cfg["runs"][run], "nickname": _nick(run)}
            for inst in ncfg["instruments"]:
                for product in ncfg["products"]:
                    render_native(master, run, inst, product, ncfg, spec, boxes)

if __name__ == "__main__":
    main()
