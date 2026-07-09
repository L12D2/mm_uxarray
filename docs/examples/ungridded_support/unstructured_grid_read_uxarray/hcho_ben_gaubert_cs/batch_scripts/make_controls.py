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

def main():
    cfg = yaml.safe_load(open(HERE / "runs.yaml"))
    master = yaml.safe_load(open(HERE / cfg.get("master", "control_master.yaml")))
    common = cfg.get("common")        
    wanted = sys.argv[1:] or list(cfg["runs"])
    for name in wanted:
        render(master, name, cfg["runs"][name], common=common)


if __name__ == "__main__":
    main()
