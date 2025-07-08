
"""
Credit the code to ana-v-espinoza
# # https://github.com/ana-v-espinoza/tropopauseCalc/blob/main/wmo.py

Description: 
# # this is the generic version listed in ISSUE 334.
# # looks like there are a couple different options that could be incorporated: 
# # https://github.com/ana-v-espinoza/tropopauseCalc/blob/main/tropCalc.py

# current issues:
Runs a tropopause calculation for each grid cell the sonde crosses -- I think. 

Consistently getting "low tropopause" (e.g. < 85 hPa/8500 Pa). This could be due to unit errors. 

#################### current error if ran as is: 
# ValueError: Shape for z_src (11, 64, 1) is not a subset of fz_src (11, 1).

A tropopause obs calculation has yet to be developed. Initial scripts can be found in beta_development. 

Author: Liam Thompson
__status__ = under development. 

"""

import numpy as np
import xarray as xr
from metpy.calc import thickness_hydrostatic, sigma_to_pressure
from metpy.units import units

# # # incorporate height / pressure coords

def wmo_tropo(obj, varmap=None, output_key="wmo_tropo_mod", lapseC=2.0, return_height=False):

    # Variable mapping
    temperature_key = varmap["temp_mod"] if varmap and "temp_mod" in varmap else "temperature_k"
    pressure_key = varmap["pres_mod"] if varmap and "pres_mod" in varmap else "pressure_model"

    pressure = obj[pressure_key] # must be in Pa
    temperature = obj[temperature_key] # must be in K

    # Identify vertical dimension
    non_vert_dims = {"time", "x", "y"}
    vert_dim = next((dim for dim in pressure.dims if dim not in non_vert_dims), None)
    if vert_dim is None:
        raise ValueError("Could not determine vertical dimension for pressure.")

    print("Available variables:", list(obj.data_vars))

    # Constants
    g = 9.80665  # m/s^2
    R = 287.05   # J/(kg·K)
    const = g / R  # 1/K/m
    pMin = 8500.0  # Pa
    pMax = 45000.0  # Pa
    dZ = 2000.0  # meters
    
    # WMO tropopause logic
    def wmo_profile(p, T):
        nLev = p.size
        nLevm = nLev - 1
        found = False

        lapse = np.zeros(nLevm)
        pHalf = np.zeros(nLevm)
        pTrop = 0.0
        
        for i in range(nLevm):
            lapse[i] = const * np.log(T[i] / T[i + 1]) / np.log(p[i] / p[i + 1])
            pHalf[i] = 0.5 * (p[i] + p[i + 1])

        for i in range(nLevm - 1):
            if lapse[i] < lapseC / 1000.0 and p[i] < pMax and not found:
                P1 = np.log(pHalf[i])
                P2 = np.log(pHalf[i + 1])
                if lapse[i] != lapse[i + 1]:
                    weight = (lapseC / 1000.0 - lapse[i]) / (lapse[i + 1] - lapse[i])
                    pTrop = np.exp(P1 + weight * (P2 - P1))
                else:
                    pTrop = pHalf[i]

                p2km = pTrop * np.exp(-dZ * const / T[i])
                lapseSum = 0
                kount = 0
                for L in range(i, nLevm):
                    if pHalf[L] > p2km:
                        lapseSum += lapse[L]
                        kount += 1
                lapseAvg = lapseSum / kount if kount > 0 else lapse[i]
                found = lapseAvg < lapseC / 1000.0
                
                if found:
                    if pTrop < pMin:
                        print(f"Warning: Tropopause pressure unusually low ({pTrop:.2f} Pa)")
                    return pTrop  # in Pa

        if return_height:
            z = thickness_hydrostatic(p[0:i], T[0:i])
            print("z in km:", z)
            return z  # returns a Quantity object with units
        return pTrop  # returns pressure in Pascals

    # Apply using xarray.apply_ufunc
    trop = xr.apply_ufunc(
        wmo_profile,
        pressure,
        temperature,
        input_core_dims=[[vert_dim], [vert_dim]],
        output_core_dims=[[]],
        vectorize=True,
        dask="parallelized",
        dask_gufunc_kwargs={"allow_rechunk": True},
        output_dtypes=[float],
    )

    # Assign to dataset
    output_dims = tuple(d for d in pressure.dims if d != vert_dim)
    # trop_da = xr.DataArray(
    #     data=trop.data,
    #     dims=output_dims,
    #     coords={dim: obj.coords[dim] for dim in output_dims if dim in obj.coords},
    #     attrs={"units": "km" if return_height else "mbar"}
    # )

    # set up the code to save output as a .csv for further debugging. 
    coords = {dim: obj.coords[dim] for dim in output_dims if dim in obj.coords}
    
    # Explicitly add lat/lon
    for coord in ["lat", "latitude", "lon", "longitude"]:
        if coord in obj.coords:
            coords[coord] = obj.coords[coord]
    
    trop_da = xr.DataArray(
        data=trop.data,
        dims=output_dims,
        coords=coords,
        attrs={"units": "m" if return_height else "Pa"}
    )

    obj[output_key] = trop_da
    
    #print(type(obj))
    #print(trop_da)
    
    return obj

    #################### current error if ran as is: 
    # ValueError: Shape for z_src (11, 64, 1) is not a subset of fz_src (11, 1).

    ######### Incorporate an option to save output as a .csv? 
    # result = wmo_tropo(ds, varmap=varmap)
    # result["wmo_tropo_mod"].to_dataframe(name="tropopause_pressure").reset_index().to_csv("tropopause_output.csv", index=False)

    