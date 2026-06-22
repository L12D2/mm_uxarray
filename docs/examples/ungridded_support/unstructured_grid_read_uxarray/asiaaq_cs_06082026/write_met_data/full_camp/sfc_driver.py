
from melodies_monet import driver
an = driver.analysis()
an.control = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/write_met_data/full_camp/sfc_pair.yaml"

# run mm
an.read_control()
an.open_models()
an.open_obs()
an.pair_data()

an.plotting()

an.save_analysis()
