from melodies_monet import driver
from melodies_monet.util import obs2obs_util

an = driver.analysis()
an.control = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/obs2obs_testing/control_tempo_l2_hcho_cesm_se.yaml"
an.read_control()
an.read_analysis()                     # loads all labels in read.paired.filenames
obs2obs_util.run(an.paired, an.control_dict["obs2obs"], default_outdir=an.output_dir)
