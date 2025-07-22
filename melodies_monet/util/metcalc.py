
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

# calc potential temperature
def ptemp(obj, varmap=None, output_key="ptemp", default_keys=None):
    if default_keys is None:
        default_keys = {
            "pressure": "pressure",
            "temperature": "temperature"
        }

    pressure_key = varmap.get("pres", default_keys["pressure"]) if varmap else default_keys["pressure"]
    temperature_key = varmap.get("temperature", default_keys["temperature"]) if varmap else default_keys["temperature"]
        
    pres = obj[pressure_key]
    temp = obj[temperature_key]
    
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

