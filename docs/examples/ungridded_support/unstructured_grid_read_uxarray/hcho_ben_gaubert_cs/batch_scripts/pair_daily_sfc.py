from melodies_monet import driver

def main():
    # run mm
    an = driver.analysis()
    an.control = "/glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/sfc_img/control_airnow_cesm_se.yaml"

    an.read_control()
    an.open_models()
    #an.open_obs()
    
    print("pairing data....", flush=True)
    #an.pair_data()
    #an.save_analysis()
    an.read_analysis()
    
    print("plotting data....", flush=True)
    an.plotting()
    # an.stats()

if __name__ == "__main__":
    main()
