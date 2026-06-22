from melodies_monet import driver
an = driver.analysis()
an.control = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/control_sfc_full_camp_cesm-se.yaml"

# run mm
an.read_control()
an.read_analysis()

# an.open_models()
# an.open_obs()

an.plotting()
an.stats()


