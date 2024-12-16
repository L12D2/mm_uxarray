from melodies_monet import region_select

import xarray as xr
import numpy as np

def test_create_custom_mask():

    np.random.seed(0)
    fake_model = 10 + 2 * np.random.randn(30,20,4)
    np.random.seed(1)
    fake_obs = 10 + 2 * np.random.randn(30, 20, 4)
    longitude = np.linspace()
    mock_data = xr.Dataset()


