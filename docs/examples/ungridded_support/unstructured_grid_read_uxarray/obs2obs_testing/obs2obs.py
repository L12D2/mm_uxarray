from melodies_monet import driver
from melodies_monet.util import obs2obs_util
import os

an = driver.analysis()
an.control = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/obs2obs_testing/control_tempo_l2_hcho_cesm_se.yaml"
an.read_control()

cfg  = an.control_dict["obs2obs"]
only = os.environ.get("OBS2OBS_ONLY") or None          # e.g. "coupling" or "em_no2"
used = obs2obs_util.labels_used(cfg, only=only)
fn   = an.control_dict["analysis"]["read"]["paired"]["filenames"]
an.control_dict["analysis"]["read"]["paired"]["filenames"] = {
    k: v for k, v in fn.items() if k in used}
print(f"obs2obs: opening {len(used)} of {len(fn)} labels", flush=True)

an.read_analysis()                     # loads all labels in read.paired.filenames

# per-label RAM footprint, so an OOM names the culprit instead of dying silent
_tot = 0.0
for _lab, _p in an.paired.items():
    try:
        _o = _p.obj
        _gb = (_o.nbytes if hasattr(_o, "nbytes")
               else _o.memory_usage(deep=True).sum()) / 1e9
    except Exception:
        _gb = float("nan")
    _tot += 0.0 if _gb != _gb else _gb          # skip NaN
    _sz = dict(_o.sizes) if hasattr(_o, "sizes") else _o.shape
    print(f"obs2obs: loaded {_lab}: {_gb:.2f} GB  sizes={_sz}", flush=True)
print(f"obs2obs: total loaded ~{_tot:.1f} GB across {len(an.paired)} labels",
      flush=True)

obs2obs_util.run(an.paired, cfg, default_outdir=an.output_dir, only=only)
