
"""Render per-run PLOTTING controls for the native-resolution TEMPO city product.
akin to make_control.py
"""
import copy
import pathlib
import os

import yaml

HERE = pathlib.Path(__file__).parent
NATIVE_ROOT = "/glade/work/lcthompson/mm_output"
# RES = 0.03
# RESTAG = ("%g" % RES).replace(".", "p")

METHOD_TAG = {"conservative": "cons", "conservative_normed": "consn",
              "radius_mean": "rmean", "nearest_s2d": "nn", "nearest_d2s": "nnd",
              "bilinear": "bilin"}
METHOD = os.environ.get("MM_NATIVE_METHOD", "conservative").strip()

MTAG = METHOD_TAG.get(METHOD, METHOD[:5])
MIN_OBS = int(os.environ.get("MM_NATIVE_MIN_OBS", "3"))

# "obsgrid": lat/lon L3 map (grp_native_* spatial_overlay plots)
# "swath":  native pixel vector, plotted standalone via plot_sat.py PLOT_SWATH

PRODUCT = os.environ.get("MM_NATIVE_PRODUCT", "obsgrid").strip().lower()
assert PRODUCT in ("obsgrid", "swath"), f"MM_NATIVE_PRODUCT={PRODUCT!r}"

INSTRUMENT = os.environ.get("MM_NATIVE_INSTRUMENT", "tempo").strip().lower()
assert INSTRUMENT in ("tempo", "tropomi"), f"MM_NATIVE_INSTRUMENT={INSTRUMENT!r}"

RES = 0.03 if INSTRUMENT == "tempo" else 0.05     # must match pair_*_native.py res
RESTAG = ("%g" % RES).replace(".", "p")

SPECIES_BY_INSTRUMENT = {
    "tempo":   {"hcho": "tempo_l2_hcho", "no2": "tempo_l2_no2"},
    "tropomi": {"no2": "tropomi_l2_no2", "hcho": "tropomi_l2_hcho",
                "co": "tropomi_l2_co"},
}

# cities paired per run (must match what you actually ran through pair_tempo_native.py)
RUN_CITIES = {
    "nonbiog": ["atl", "dfw", "la", "den"],
    "biog":    ["atl", "dfw", "la", "den"],
    "grapes":  ["atl", "dfw", "la", "den"],
    "mxcat":   ["atl", "mex"],
}
# short city id -> box name used in runs.yaml common.add_domains
CITY_BOX = {"atl": "ATL_metro", "mex": "MEX_metro", "la": "LA_metro",
            "den": "DEN_metro", "dfw": "DFW_metro", "seus": "SE_US"}

SPECIES = SPECIES_BY_INSTRUMENT[INSTRUMENT]

def main():
    runs_cfg = yaml.safe_load(open(HERE / "runs.yaml"))
    boxes = (runs_cfg.get("common", {}).get("add_domains") or {})
    skel = yaml.safe_load(open(HERE / "control_tempo_native.yaml"))

    for run, cities in RUN_CITIES.items():
        rspec = runs_cfg["runs"][run]
        paired_dir = f"{NATIVE_ROOT}/{run}_{INSTRUMENT}_native"
        plot_dir = rspec["plot_dir"].rstrip("/") + (
            "_native" if INSTRUMENT == "tempo" else f"_{INSTRUMENT}_native")

        cd = copy.deepcopy(skel)
        an = cd["analysis"]
        an["output_dir"] = plot_dir
        an["output_dir_save"] = paired_dir
        an["output_dir_read"] = paired_dir
        an["montage"] = {
            "plot_dir": plot_dir, "outdir": plot_dir + "/montages",
            "group_by": "type", "search": "*.png", "cols": 6, "thumb_width": 420,
        }
        an.pop("save", None)  # plotting only

        cd["model"]["cam-chem-se"]["files"] = rspec["model_files"]
        cd["model"]["cam-chem-se"]["scrip_file"] = rspec["scrip_file"]

        filenames, plots = {}, {}
        for city in cities:
            box = CITY_BOX[city]
            if box not in boxes:
                print(f"  WARNING: {run}/{city}: box {box} not in runs.yaml; skipping.")
                continue
            for sp, obs_name in SPECIES.items():
                label = f"{city}_{obs_name}_cam-chem-se_{PRODUCT}"
                filenames[label] = [
                    f"{paired_dir}/{city}_{RESTAG}_{MTAG}_*_{obs_name}_cam-chem-se_{PRODUCT}.nc4"
                ]

                # shared, generalized plot-group skeleton (no hardcoding per city)
                common = {
                    "domain_type": ["auto-region:box"],
                    "domain_name": [box],
                    "domain_info": {box: boxes[box]},
                    "data": [label],
                }
                if PRODUCT == "swath":
                    # native pixel vector -> generalized spatial_swath type,
                    # dispatched through an.plotting() (driver _plot_spatial_swath).
                    plots[f"grp_swath_{city}_{sp}"] = {
                        "type": "spatial_swath",
                        "fig_kwargs": {"figsize": [22, 6], "states": True},
                        "text_kwargs": {"fontsize": 14},
                        "data_proc": {"render": "auto"},   # bin_deg auto by density
                        **common,
                    }
                else:  # obsgrid -> gridded spatial_overlay (unchanged)
                    plots[f"grp_native_{city}_{sp}"] = {
                        "type": "spatial_overlay",
                        "fig_kwargs": {"figsize": [18, 5], "cbar_orientation": "horizontal",
                                       "cbar_kwargs": {"shrink": 0.6},
                                       "states": True, "counties": True},
                        "text_kwargs": {"fontsize": 16},
                        "data_proc": {"time_reduction": "mean", "daily_first": True,
                                      "common_mask": True, "min_obs": MIN_OBS,
                                      "set_axis": False, "rem_obs_nan": True},
                        **common,
                    }
                    
        an["read"] = {"paired": {"method": "netcdf", "filenames": filenames}}
        cd["plots"] = plots
        cd.pop("stats", None)

        _psfx = "" if PRODUCT == "obsgrid" else f"_{PRODUCT}"
        _itag = "" if INSTRUMENT == "tempo" else f"{INSTRUMENT}_"
        out = HERE / f"control_{run}_{_itag}native_{MTAG}{_psfx}.yaml"
        with open(out, "w") as f:
            f.write(f"# GENERATED by make_native_controls.py "
                    f"(instrument: {INSTRUMENT}, run: {run}, product: {PRODUCT}).\n"
                    "# Do not edit -- edit runs.yaml / control_tempo_native.yaml and regenerate.\n")
            yaml.safe_dump(cd, f, sort_keys=False, width=120)
        print(f"wrote {out}  ({len(plots)} groups, {len(filenames)} labels)")


if __name__ == "__main__":
    main()