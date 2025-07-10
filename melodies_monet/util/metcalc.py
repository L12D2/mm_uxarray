
"""
This file calculates a number of meteorological variables. For hydrometeorological variables, 
variables including dewpoint and relative humidity are not universally provided by the obs. 
For models, wind speed is usually not directly provided. 

This script calculates all of them in one place. 

Users can calculate variables directly from the observations by using the extra_calc dictionary. Potential 
temperature is the only variable with this supported capability. 

See: control_aircraft_looping_AEROMMA_UFSAQM.yaml

Author: Liam Thompson
"""

# will need to make this an optional dependency if we proceed in using this. 
import metpy

# import specific metpy libraries needed for each calculation
from metpy.calc import dewpoint_from_specific_humidity
from metpy.calc import relative_humidity_from_specific_humidity
from metpy.calc import potential_temperature
from metpy.calc import wind_speed

# if future iterations want to calc dewpoint/relh on the observations, this can be done with other metpy
# libraries. 

# addtl libraries to make the world go round
from metpy.units import units
import metpy.constants as mconst
import numpy as np
import pandas as pd
import xarray as xr

# calc dewpoint 
def dewpoint(obj, varmap = None, output_key = "dewpoint"):
    # grab variable names from the yaml 
    pressure_key = varmap['pressure'] if varmap and 'pressure' in varmap else 'surfpres_pa'
    specific_hum_key = varmap['specific_hum'] if varmap and 'specific_hum' in varmap else 'specific_hum'

    pressure = obj[pressure_key]
    specific_hum = obj[specific_hum_key]

    dpt = (metpy.calc.dewpoint_from_specific_humidity(
        pressure * units.Pa,
        specific_hum * units("kg/kg")
    )).metpy.convert_units("K")

    dpt_np = dpt.astype("float64").values

    # DEBUG
    # make sure this is generic enough
    # print("pressure dims:", pressure.dims)
    # print("specific_hum dims:", specific_hum.dims)
    # print("dpt_np shape:", dpt_np.shape)

    # fix is needed in order to work with vert plots. Check dimensions before proceeding. 
    if hasattr(obj, "coords") and hasattr(obj, "dims"):
        # Use pressure dims if number of dims matches dpt_np
        if len(pressure.dims) == dpt_np.ndim:
            dims_to_use = pressure.dims
        # use 4d case for vert plots
        elif len(specific_hum.dims) == dpt_np.ndim:
            dims_to_use = specific_hum.dims
        else:
            raise ValueError(
                f"No matching dims for output: pressure.dims={pressure.dims}, "
                f"specific_hum.dims={specific_hum.dims}, dpt_np.shape={dpt_np.shape}"
            )

        obj[output_key] = (dims_to_use, dpt_np)
        obj[output_key].attrs["units"] = "K"
        return obj
        
    elif isinstance(obj, dict):
        obj[output_key] = dpt_np
        return obj
    else:
        return dpt_np
        
# calc relative humidity
def relh(obj, varmap=None, output_key="rel_hum"):
    # grab variable names from the yaml or fall back to defaults
    pressure_key = varmap["pressure"] if varmap and "pressure" in varmap else "surfpres_pa"
    specific_hum_key = varmap["specific_hum"] if varmap and "specific_hum" in varmap else "specific_hum"
    temperature_key = varmap["temperature"] if varmap and "temperature" in varmap else "temperature_k"

    pressure = obj[pressure_key]
    specific_hum = obj[specific_hum_key]
    temperature = obj[temperature_key]

    rlh = (metpy.calc.relative_humidity_from_specific_humidity(
        pressure * units.Pa,
        temperature * units.kelvin, 
        specific_hum * units("kg/kg")
    )).metpy.convert_units("%")

    rlh_np = rlh.astype("float64").values

    # Set any value above 100% to NaN
    # I think this happens because the divisions by 0 
    # (e.g. very tinyyyyyy numbers once you get high enough in the atmos) 
    # are handled  weirdly in the metpy library. 
    rlh_np = np.where(rlh_np > 100, np.nan, rlh_np)

    # DEBUG
    # print("pressure dims:", pressure.dims)
    # print("specific_hum dims:", specific_hum.dims)
    # print("temperature dims:", temperature.dims, temperature.shape)
    # print("rlh_np shape:", rlh_np.shape)

    # fix is needed in order to work with vert plots. Check dimensions before proceeding. 
    if hasattr(obj, "coords") and hasattr(obj, "dims"): 
        dims_to_use = list(specific_hum.dims)
        if rlh_np.shape != specific_hum.shape:
            try:
                axes_needed = [rlh_np.shape.index(s) for s in specific_hum.shape]
                rlh_np = np.transpose(rlh_np, axes_needed)
                #print("Transposed rlh_np shape:", rlh_np.shape)
            except Exception as e:
                #print("Error transposing rlh_np:", e)
                raise ValueError(
                    f"Could not align output shape {rlh_np.shape} to input shape {specific_hum.shape}; "
                    "please check dimension alignment."
                )
        obj[output_key] = (dims_to_use, rlh_np)
        obj[output_key].attrs["units"] = "%"
        return obj
    elif isinstance(obj, dict):
        obj[output_key] = rlh_np
        return obj
    else:
        return rlh_np
        
# calc windspeed
def wspd(obj, varmap = None, output_key = "windspeed"):
    # grab variable names from the yaml 
    u_key = varmap["u_comp"] 
    v_key = varmap["v_comp"] 

    u = obj[u_key]
    v = obj[v_key]

    wspd = (metpy.calc.wind_speed(
        u * units("m/s"),
        v * units("m/s")
        )).metpy.convert_units("m/s")

    wspd_np = wspd.astype("float64").values

    # DEBUG
    # print("u dims:", u.dims)
    # print("v dims:", v.dims)
    # print("wspd_np shape:", wspd_np.shape)

    # fix is needed in order to work with vert plots. Check dimensions before proceeding. 
    if hasattr(obj, "coords") and hasattr(obj, "dims"): 
        dims_to_use = list(u.dims)
        if wspd_np.shape != u.shape:
            try:
                axes_needed = [wspd_np.shape.index(s) for s in u.shape]
                wspd_np = np.transpose(wspd_np, axes_needed)
                #print("Transposed wspd_np shape:", wspd_np.shape)
            except Exception as e:
                print("Error transposing wspd_np:", e)
                raise ValueError(
                    f"Could not align output shape {wspd_np.shape} to input shape {u.shape}; "
                    "please check dimension alignment."
                )
        obj[output_key] = (dims_to_use, wspd_np)
        obj[output_key].attrs["units"] = "m/s"
        return obj
    elif isinstance(obj, dict):
        obj[output_key] = wspd_np
        return obj
    else:
        return wspd_np

# calc wind direction
def wdir(obj, varmap = None, output_key = "winddir"):
    # grab variable names from the yaml 
    u_key = varmap["u_comp"] 
    v_key = varmap["v_comp"] 
    
    u = obj[u_key]
    v = obj[v_key]

    # metpy version of this is throwing in weird dimensions so calc by hand
    # wdr_rad = np.arctan2(-v, -u)
    wdr_rad = np.arctan2(u, v)
    #wdr_deg = np.degrees(wdr_rad)
    wdr_deg = (np.degrees(wdr_rad) + 180)
    winddir = wdr_deg % 360
    #print(winddir)
    winddir_np = winddir.astype("float64").values

    if hasattr(obj, "coords") and hasattr(obj, "dims"):  
        dims_to_use = list(u.dims)
        if winddir_np.shape != u.shape:
            try:
                axes_needed = [winddir_np.shape.index(s) for s in u.shape]
                winddir_np = np.transpose(winddir_np, axes_needed)
                #print("Transposed winddir_np shape:", winddir_np.shape)
            except Exception as e:
                print("Error transposing winddir_np:", e)
                raise ValueError(
                    f"Could not align output shape {winddir_np.shape} to input shape {u.shape}; "
                    "please check dimension alignment."
                )
        obj[output_key] = (dims_to_use, winddir_np)
        obj[output_key].attrs["units"] = "degrees"
        return obj
    elif isinstance(obj, dict):
        obj[output_key] = winddir_np
        return obj
    else:
        return winddir_np

# calc potential temperature
def ptemp(obj, varmap=None, output_key="ptemp", default_keys=None):
    if default_keys is None:
        default_keys = {
            "pressure": "pressure",
            "temperature": "temperature"
        }

    pressure_key = varmap.get("pres", default_keys["pressure"]) if varmap else default_keys["pressure"]
    temperature_key = varmap.get("temp", default_keys["temperature"]) if varmap else default_keys["temperature"]
        
    pres = obj[pressure_key]
    temp = obj[temperature_key]

    # # # the unit conversion for celsius for the obs may not be occuring till later
    # if temp.max() < 200:
    #     print(f"Detected {output_key} temperature likely in Celsius. Converting to Kelvin.")
    #     temp = temp + 273.15
    #     print(temp)
    # else:
    #     print("continue")
    #     #print(f"{output_key} temp values: {temp.isel(z=slice(0, 5)).values}")
    
    ptmp = (metpy.calc.potential_temperature(
        pres * units.Pa,
        temp * units("K")
    )).metpy.convert_units("K")

    ptmp_np = ptmp.astype("float64").values
    
    if hasattr(obj, "coords") and hasattr(obj, "dims"):
        if len(pres.dims) == ptmp_np.ndim:
            dims_to_use = pres.dims
        elif len(temp.dims) == ptmp_np.ndim:
            dims_to_use = temp.dims
        else:
            raise ValueError(
                f"No matching dims for output: pressure.dims={pres.dims}, "
                f"temp.dims={temp.dims}, ptmp_np.shape={ptmp_np.shape}"
            )

        obj[output_key] = (dims_to_use, ptmp_np)
        obj[output_key].attrs["units"] = "K"
        return obj

    elif isinstance(obj, dict):
        obj[output_key] = ptmp_np
        return obj
    else:
        return ptmp_np

# calc bulk richardson number:
# https://github.com/ThomasRieutord/MetPy/blob/add_boundarylayer_module/src/metpy/calc/boundarylayer.py
def bulk_richardson_number(
    height,
    potential_temperature,
    u,
    v,
    idxfoot: int = 0,
    ustar=0 * units.meter_per_second,
):
    r"""Calculate the bulk Richardson number.

    See [VH96], eq. (3):

    .. math::   Ri = (g/\theta) * \frac{(\Delta z)(\Delta \theta)}
             {\left(\Delta u)^2 + (\Delta v)^2 + b(u_*)^2}

    Parameters
    ----------
    height : `pint.Quantity`
        Altitude (metres above ground) of the points in the profile
    potential_temperature : `pint.Quantity`
        Potential temperature profile
    u : `pint.Quantity`
        Zonal wind profile
    v : `pint.Quantity`
        Meridional wind profile
    idxfoot : int, optional
        The index of the foot point (first trusted measure), defaults to 0.

    Returns
    -------
    `pint.Quantity`
        Bulk Richardson number profile
    """
    if idxfoot == 0:
        # Force the ground level to have null wind
        Du = u
        Dv = v
    else:
        Du = u - u[idxfoot]
        Dv = v - v[idxfoot]
    
    Dtheta = potential_temperature - potential_temperature[idxfoot]
    Dz = height - height[idxfoot]

    idx0 = Du**2 + Dv**2 + ustar**2 == 0
    if idx0.sum() > 0:
        bRi = np.ones_like(Dtheta) * np.nan * units.dimensionless
        bRi[~idx0] = (
            (mconst.g / potential_temperature[~idx0])
            * (Dtheta[~idx0] * Dz[~idx0])
            / (Du[~idx0] ** 2 + Dv[~idx0] ** 2 + ustar**2)
        )
    else:
        bRi = (
            (mconst.g / potential_temperature)
            * (Dtheta * Dz)
            / (Du**2 + Dv**2 + ustar**2)
        )

    return bRi

# calc blh from richardson bulk 
# https://github.com/ThomasRieutord/MetPy/blob/add_boundarylayer_module/src/metpy/calc/boundarylayer.py
def blh_from_richardson_bulk(
    height,
    potential_temperature,
    u,
    v,
    smoothingspan: int = 10,
    idxfoot: int = 0,
    bri_threshold=0.25 * units.dimensionless,
    ustar=0.1 * units.meter_per_second,
):
    """Calculate atmospheric boundary layer height with the method of
    bulk Richardson number.

    It is the height where the bulk Richardson number exceeds a given threshold.
    Well indicated for unstable boundary layers. See [VH96, Sei00, Col14, Guo16].

    Parameters
    ----------
    height : `pint.Quantity`
        Altitude (metres above ground) of the points in the profile
    potential_temperature : `pint.Quantity`
        Potential temperature profile
    u : `pint.Quantity`
        Zonal wind profile
    v : `pint.Quantity`
        Meridional wind profile
    smoothingspan : int, optional
        The amount of smoothing (number of points in moving average)
    idxfoot : int, optional
        The index of the foot point (first trusted measure), defaults to 0.
    bri_threshold : `pint.Quantity`, optional
        Threshold to exceed to get boundary layer top. Defaults to 0.25
    ustar : `pint.Quantity`, optional
        Additional friction term in [VH96]. Defaluts to 0.

    Returns
    -------
    blh : `pint.Quantity`
        Boundary layer height estimation
    """

    # Apply rolling mean smoothing via pandas
    temp_smooth = pd.Series(potential_temperature.magnitude).rolling(
        window=smoothingspan, center=True, min_periods=1
    ).mean().values * units.kelvin

    u_smooth = pd.Series(u.magnitude).rolling(
        window=smoothingspan, center=True, min_periods=1
    ).mean().values * u.units

    v_smooth = pd.Series(v.magnitude).rolling(
        window=smoothingspan, center=True, min_periods=1
    ).mean().values * v.units

    # Calculate Richardson number profile
    bRi = bulk_richardson_number(
        height, temp_smooth, u_smooth, v_smooth,
        idxfoot=idxfoot, ustar=ustar
    )

    # Mask NaNs
    mask = ~np.isnan(bRi)
    height = height[mask]
    bRi = bRi[mask]

    # Find first level above threshold
    if np.any(bRi > bri_threshold):
        iblh = np.where(bRi > bri_threshold)[0][0]
        blh = height[iblh]
    else:
        blh = np.nan * units.meter

    return blh
    
# calc / estimate boundary layer
# def blayer(df):
#     """
#     Description: This function estimates the boundary layer height using three different methods. 
#            1) Potential temperature inversion 
#            2) Concentration gradient (specific humidity and mixing ratio)
#            3) Richardson numbers 

#            Note: These are estimated heights. Users should evaluate these calculations against the 
#            vertprofile sonde. 
           
#     Source and credit for inspiration and guidance on building this code:
#     https://github.com/ThomasRieutord/MetPy/blob/add_boundarylayer_module/src/metpy/calc/boundarylayer.py


#     Addtl information that may be of use: 
#     https://www.weather.gov/media/zhu/ZHU_Training_Page/clouds/planetary_boundary_layer/L1-PBL.pdf
#     """
    
#     # pull appropriate cols from df
#     ptemp_obs = df["ptemp_obs"].values * units.kelvin
#     ptemp_mod = df["ptemp_mod"].values * units.kelvin
#     height = df["ghght_obs"].values * units.meter
#     spfh_kgkg = df["spfh"].values * units("kg/kg")
#     mxng_ratio_kgkg = df["w"].values * units("kg/kg")
    
#     # u_comp = extra_calc.get('ptemp_mod', {}).get('u_comp', None)
#     # v_comp = extra_calc.get('ptemp_mod', {}).get('v_comp', None)
#     u_comp = df["ugrd"].values * units.meter / units.second
#     v_comp = df["vgrd"].values * units.meter / units.second
    
#     # skip the lowest 100 m. Noisy.
#     idxfoot = np.argmax(height > 100 * units.meter)
#     height = height[idxfoot: ] 
#     ptemp_obs = ptemp_obs[idxfoot: ] 
#     ptemp_mod = ptemp_mod[idxfoot: ]
#     u_comp = u_comp[idxfoot: ]
#     v_comp = v_comp[idxfoot: ]
#     spfh_kgkg = spfh_kgkg[idxfoot: ]
#     mxng_ratio_kgkg = mxng_ratio_kgkg[idxfoot: ]
    
#     # sort the trim
#     sort_idx = np.argsort(height)
#     height = height[sort_idx]
#     ptemp_obs = ptemp_obs[sort_idx]
#     ptemp_mod = ptemp_mod[sort_idx]
#     u_comp = u_comp[sort_idx]
#     v_comp = v_comp[sort_idx]
#     spfh_kgkg = spfh_kgkg[sort_idx]
#     mxng_ratio_kgkg = mxng_ratio_kgkg[sort_idx]
    
#     # smooth the ptemp to help reduce noise. Pandas is NaN aware. 
#     size = 3
#     ptemp_smooth_obs = pd.Series(ptemp_obs.magnitude).rolling(window=size, 
#                                                               center=True, min_periods=1).mean().values * units.kelvin
#     ptemp_smooth_mod = pd.Series(ptemp_mod.magnitude).rolling(window=size, 
#                                                               center=True, min_periods=1).mean().values * units.kelvin
    
#     dtheta_dz_obs = metpy.calc.first_derivative(ptemp_smooth_obs, x=height)
#     dtheta_dz_mod = metpy.calc.first_derivative(ptemp_smooth_mod, x=height)

#     # debug 
#     # print(dtheta_dz_obs)
#     # print(dtheta_dz_mod)

#     # set a minimum gradient and altitude in which you would normally expect to find the BL
#     min_grad = 0.005 * units("K/m")
#     min_alt = 750 * units.meter

#     # boundary layer can be found when d theta/dz > 0 .
#     # this can be sensitive due to random noise. So, a threshold is needed. 
#     # update to make a list of thresholds to test the BL at for the temp inversion method. 
#     def detect_blh(dtheta_dz, height):
#         valid = (dtheta_dz > min_grad) & (height > min_alt)
#         if np.any(valid):
#             idx = np.where(valid)[0][0]
#             return height[1:][idx]
        
#         if dtheta_dz < min_grad or height < min_alt:
#             print("Based on a minimum altitude of 500 m and minimum gradient of detection, the boundary layer approximation should be used with caution.")
#             return 

#     # # Detect BLH
#     blh_obs = detect_blh(dtheta_dz_obs, height)
#     blh_mod = detect_blh(dtheta_dz_mod, height)
    
#     #print("Observed BLH from potential temperature inversion:", blh_obs)
#     #print("Modeled BLH from potential temperature inversion:",blh_mod)

#     # concentration method. Output for specific humidity (from the model) and mixing ratio (obs)
#     size = 3
#     mxng_smooth_obs = pd.Series(mxng_ratio_kgkg.magnitude).rolling(window=size, 
#                                                               center=True, min_periods=1).mean().values * units("kg/kg")
#     spfh_smooth_mod = pd.Series(spfh_kgkg.magnitude).rolling(window=size, 
#                                                               center=True, min_periods=1).mean().values * units("kg/kg")
    
#     dcdz_spfh = metpy.calc.first_derivative(spfh_smooth_mod, x=height)
#     dcdz_mxng = metpy.calc.first_derivative(mxng_smooth_obs, x=height)
    
#     spfh_blh = np.argmin(dcdz_spfh)
#     mxng_blh = np.argmin(dcdz_mxng)

#     # Specific humidity BLH
#     if height[spfh_blh] > min_alt:
#         spfh_blh_val = height[spfh_blh]
#     else:
#         for i in range(spfh_blh + 1, len(height)):
#             if height[i] > min_alt:
#                 spfh_blh_val = height[i]
#                 break
#         else:
#             spfh_blh_val = height[-1] 

#     # Mixing ratio BLH
#     if height[mxng_blh] > min_alt:
#         mxng_blh_val = height[mxng_blh]
#     else:
#         for i in range(mxng_blh + 1, len(height)):
#             if height[i] > min_alt:
#                 mxng_blh_val = height[i]
#                 break
#         else:
#             mxng_blh_val = height[-1]  

#     # Debug
#     # print("spfh",spfh_kgkg )
#     # print(len(spfh_kgkg ))
#     # print("mix",mxng_ratio_kgkg)
#     # print(len(mxng_ratio_kgkg))
#     # print(dcdz_spfh)
#     # print(dcdz_mxng)

#     #print("Observed mixing ratio BLH:", height[mxng_blh])
#     #print("Modeled specific humidity BLH:", height[spfh_blh])

#     # richardson number method 
#     # calculate the richardson number 
#     r_num_obs = bulk_richardson_number(height = height, potential_temperature = ptemp_obs, 
#                                        u = u_comp, v = v_comp, 
#                                        idxfoot = 0, ustar=0 * units.meter_per_second)
    
#     r_num_mod = bulk_richardson_number(height = height, potential_temperature = ptemp_mod, 
#                                        u = u_comp, v = v_comp, 
#                                        idxfoot = 0, ustar=0 * units.meter_per_second)

#     # Calc blh from rich number 
#     r_blh_obs = blh_from_richardson_bulk(height = height, potential_temperature = ptemp_obs,
#                                         u = u_comp, v = v_comp,
#                                         smoothingspan = 10,
#                                         idxfoot = 0, 
#                                         bri_threshold=0.25 * units.dimensionless,
#                                         ustar=0.1 * units.meter_per_second)

#     r_blh_mod = blh_from_richardson_bulk(height = height, potential_temperature = ptemp_mod,
#                                         u = u_comp, v = v_comp,
#                                         smoothingspan = 10,
#                                         idxfoot = 0, 
#                                         bri_threshold=0.25 * units.dimensionless,
#                                         ustar=0.1 * units.meter_per_second)

#     # richardson obs BLH
#     def find_blh_fallback(blh, height_array, min_alt):
#         if blh is None:
#             print("Warning: BLH not found, defaulting to highest available height.")
#             return height_array[-1]
    
#         # If value already valid
#         if blh > min_alt:
#             return blh
    
#         # Get index of closest height value
#         idx = np.argmin(np.abs(height_array - blh))
    
#         # Scan forward
#         for i in range(idx + 1, len(height_array)):
#             if height_array[i] > min_alt:
#                 return height_array[i]
    
#         # Fallback if all fail
#         return height_array[-1]
        
#     r_blh_obs = find_blh_fallback(r_blh_obs, height, min_alt)
#     r_blh_mod = find_blh_fallback(r_blh_mod, height, min_alt)

#     #print(r_blh_mod_val)
#     #print(r_blh_obs_val)
#     #print("Observed richardson # BLH:", r_blh_obs)
#     #print("Modeled richardson # BLH:", r_blh_mod)

#     results = {
#     "Observed BLH from potential temperature inversion": blh_obs,
#     "Modeled BLH from potential temperature inversion": blh_mod,
#     "Observed mixing ratio BLH": height[mxng_blh],
#     "Modeled specific humidity BLH": height[spfh_blh],
#     "Observed richardson # BLH": r_blh_obs,
#     "Modeled richardson # BLH": r_blh_mod,
#     }
    
#     return results

# boundary layer with a range of heights. 
def blayer(df):
    """
    Description: This function estimates the boundary layer height using three different methods,
                 now extended to provide estimates within specified height ranges.

                 1) Potential temperature inversion
                 2) Concentration gradient (specific humidity and mixing ratio)
                 3) Richardson numbers

                 The boundary layer height in this version is tested between a number of intervals. 

    Note: These are estimated heights. Users should evaluate these calculations against the
          vertprofile sonde.

    Source and credit for inspiration and guidance on building this code:
    https://github.com/ThomasRieutord/MetPy/blob/add_boundarylayer_module/src/metpy/calc/boundarylayer.py

    Addtl information that may be of use:
    https://www.weather.gov/media/zhu/ZHU_Training_Page/clouds/planetary_boundary_layer/L1-PBL.pdf
    """

    # pull appropriate cols from df
    ptemp_obs = df["ptemp_obs"].values * units.kelvin
    ptemp_mod = df["ptemp_mod"].values * units.kelvin
    height = df["ghght_obs"].values * units.meter
    spfh_kgkg = df["spfh"].values * units("kg/kg")
    mxng_ratio_kgkg = df["w"].values * units("kg/kg")

    u_comp = df["ugrd"].values * units.meter / units.second
    v_comp = df["vgrd"].values * units.meter / units.second

    # skip the lowest 100 m. Noisy.
    idxfoot = np.argmax(height > 100 * units.meter)
    height_trimmed = height[idxfoot:]
    ptemp_obs_trimmed = ptemp_obs[idxfoot:]
    ptemp_mod_trimmed = ptemp_mod[idxfoot:]
    u_comp_trimmed = u_comp[idxfoot:]
    v_comp_trimmed = v_comp[idxfoot:]
    spfh_kgkg_trimmed = spfh_kgkg[idxfoot:]
    mxng_ratio_kgkg_trimmed = mxng_ratio_kgkg[idxfoot:]

    # sort the trim
    sort_idx = np.argsort(height_trimmed)
    height_trimmed = height_trimmed[sort_idx]
    ptemp_obs_trimmed = ptemp_obs_trimmed[sort_idx]
    ptemp_mod_trimmed = ptemp_mod_trimmed[sort_idx]
    u_comp_trimmed = u_comp_trimmed[sort_idx]
    v_comp_trimmed = v_comp_trimmed[sort_idx]
    spfh_kgkg_trimmed = spfh_kgkg_trimmed[sort_idx]
    mxng_ratio_kgkg_trimmed = mxng_ratio_kgkg_trimmed[sort_idx]

    # provide ranges the BLH could be between
    height_ranges = [
        (500, 1000),
        (1000, 1500),
        (1500, 2000),
        (2000, 2500),
        (2500, 3000),
    ]

    # height_ranges = [
    #     (100, 400),      
    #     (400, 800),      
    #     (800, 1600),     
    #     (1600, 2500),    
    #     (2500, float('inf'))
    # ]
    
    results = {}

    for h_min, h_max in height_ranges:
        range_str = f"{h_min}-{h_max}m"
        results[range_str] = {}

        # Filter data for the current height range
        in_range_indices = (height_trimmed >= h_min * units.meter) & \
                           (height_trimmed < h_max * units.meter)

        if not np.any(in_range_indices):
            results[range_str] = {
                "Observed BLH from potential temperature inversion": None,
                "Modeled BLH from potential temperature inversion": None,
                "Observed mixing ratio BLH": None,
                "Modeled specific humidity BLH": None,
                "Observed richardson # BLH": None,
                "Modeled richardson # BLH": None,
            }
            continue

        height_in_range = height_trimmed[in_range_indices]
        ptemp_obs_in_range = ptemp_obs_trimmed[in_range_indices]
        ptemp_mod_in_range = ptemp_mod_trimmed[in_range_indices]
        u_comp_in_range = u_comp_trimmed[in_range_indices]
        v_comp_in_range = v_comp_trimmed[in_range_indices]
        spfh_kgkg_in_range = spfh_kgkg_trimmed[in_range_indices]
        mxng_ratio_kgkg_in_range = mxng_ratio_kgkg_trimmed[in_range_indices]

        # Potential temperature inversion method
        size = 3
        ptemp_smooth_obs = pd.Series(ptemp_obs_in_range.magnitude).rolling(
            window=size, center=True, min_periods=1
        ).mean().values * units.kelvin
        ptemp_smooth_mod = pd.Series(ptemp_mod_in_range.magnitude).rolling(
            window=size, center=True, min_periods=1
        ).mean().values * units.kelvin

        dtheta_dz_obs = metpy.calc.first_derivative(ptemp_smooth_obs, x=height_in_range)
        dtheta_dz_mod = metpy.calc.first_derivative(ptemp_smooth_mod, x=height_in_range)

        min_grad = 0.005 * units("K/m")

        def detect_blh_in_range(dtheta_dz, height_data, current_h_min):
            # We already filtered by h_min, so just check the gradient
            valid = (dtheta_dz > min_grad)
            if np.any(valid):
                idx = np.where(valid)[0][0]
                return height_data[idx]
            return None

        blh_obs_ptemp = detect_blh_in_range(dtheta_dz_obs, height_in_range, h_min * units.meter)
        blh_mod_ptemp = detect_blh_in_range(dtheta_dz_mod, height_in_range, h_min * units.meter)

        # Concentration method
        mxng_smooth_obs = pd.Series(mxng_ratio_kgkg_in_range.magnitude).rolling(
            window=size, center=True, min_periods=1
        ).mean().values * units("kg/kg")
        spfh_smooth_mod = pd.Series(spfh_kgkg_in_range.magnitude).rolling(
            window=size, center=True, min_periods=1
        ).mean().values * units("kg/kg")

        dcdz_spfh = metpy.calc.first_derivative(spfh_smooth_mod, x=height_in_range)
        dcdz_mxng = metpy.calc.first_derivative(mxng_smooth_obs, x=height_in_range)

        spfh_blh_idx = np.argmin(dcdz_spfh)
        mxng_blh_idx = np.argmin(dcdz_mxng)

        spfh_blh_val = height_in_range[spfh_blh_idx]
        mxng_blh_val = height_in_range[mxng_blh_idx]

        # Richardson number method
        # Ensure there's enough data for Richardson number calculation
        if len(height_in_range) > 1:
            r_blh_obs = blh_from_richardson_bulk(height=height_in_range, potential_temperature=ptemp_obs_in_range,
                                                u=u_comp_in_range, v=v_comp_in_range, smoothingspan=10,
                                                idxfoot=0, bri_threshold=0.25 * units.dimensionless,
                                                ustar=0.1 * units.meter_per_second)

            r_blh_mod = blh_from_richardson_bulk(height=height_in_range, potential_temperature=ptemp_mod_in_range,
                                                u=u_comp_in_range, v=v_comp_in_range, smoothingspan=10,
                                                idxfoot=0, bri_threshold=0.25 * units.dimensionless,
                                                ustar=0.1 * units.meter_per_second)
        else:
            r_blh_obs = None
            r_blh_mod = None

        # Store results for the current range
        results[range_str] = {
            "Observed BLH from potential temperature inversion": blh_obs_ptemp,
            "Modeled BLH from potential temperature inversion": blh_mod_ptemp,
            "Observed mixing ratio BLH": mxng_blh_val,
            "Modeled specific humidity BLH": spfh_blh_val,
            "Observed richardson # BLH": r_blh_obs,
            "Modeled richardson # BLH": r_blh_mod,
        }

    return results
    
    
