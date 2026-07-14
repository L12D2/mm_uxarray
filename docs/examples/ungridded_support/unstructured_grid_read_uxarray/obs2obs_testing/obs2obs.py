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
obs2obs_util.run(an.paired, cfg, default_outdir=an.output_dir, only=only)
