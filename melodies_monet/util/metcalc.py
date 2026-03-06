
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
# if future iterations want to calc dewpoint/relh on the observations, this can be done with other metpy
# libraries. 
# addtl libraries to make the world go round
from metpy.units import units

# calc dewpoint 
def dewpoint(obj, varmap = None, output_key = "dewpoint"):
    """Calculates dewpoint in K
    
    Parameters
    ----------
    obj : xarray dataset
        Model/obs data
    varmap : dictionary
        Dictionary defining the data column names to use
        Please provide variable names for "pres_calc" in Pa and "specific_hum" in kg/kg
    output_key : string
        String to name the new dewpoint ouput
    Returns
    -------
    xarray dataset
        Xarray dataset with applied calculation as a new data array
        
    """
    
    # grab variable names from the yaml 
    pressure_key = varmap['pres_calc'] if varmap and 'pres_calc' in varmap else 'surfpres_pa'
    specific_hum_key = varmap['specific_hum'] if varmap and 'specific_hum' in varmap else 'specific_hum'

    pressure = obj[pressure_key]
    specific_hum = obj[specific_hum_key]

    dpt = (metpy.calc.dewpoint_from_specific_humidity(
        pressure * units.Pa,
        specific_hum * units("kg/kg")
    )).metpy.convert_units("K").metpy.dequantify()
    #Fancy dequantify method puts back in standard xarray format

    obj[output_key] = dpt
        
    return obj

        
# calc relative humidity
def relh(obj, varmap=None, output_key="rel_hum"):
    """Calculates Relative Humidity in %
    
    Parameters
    ----------
    obj : xarray dataset
        Model/obs data
    varmap : dictionary
        Dictionary defining the data column names to use
        Please provide variable names for "pres_calc" in Pa, "temp_calc" in K, "specific_hum" in kg/kg
    output_key : string
        String to name the new relative humidity ouput
    Returns
    -------
    xarray dataset
        Xarray dataset with applied calculation as a new data array
        
    """
    
    # grab variable names from the yaml or fall back to defaults
    pressure_key = varmap['pres_calc'] if varmap and 'pres_calc' in varmap else 'surfpres_pa'
    specific_hum_key = varmap["specific_hum"] if varmap and "specific_hum" in varmap else "specific_hum"
    temperature_key = varmap["temp_calc"] if varmap and "temp_calc" in varmap else "temperature_k"

    pressure = obj[pressure_key]
    specific_hum = obj[specific_hum_key]
    temperature = obj[temperature_key]

    rlh = (metpy.calc.relative_humidity_from_specific_humidity(
        pressure * units.Pa,
        temperature * units.kelvin, 
        specific_hum * units("kg/kg"),
        phase = 'auto' #Needed for cold temperatures in stratosphere
    )).metpy.convert_units("%").metpy.dequantify()
    #Fancy dequantify method puts back in standard xarray format

    # Set any value above 100% to NaN
    # I think this happens because the divisions by 0 
    # (e.g. very tinyyyyyy numbers once you get high enough in the atmos) 
    # are handled  weirdly in the metpy library. 
    rlh = rlh.where(rlh <= 100)

    obj[output_key] = rlh
        
    return obj
        
# calc windspeed
def wspd(obj, varmap = None, output_key = "windspeed"):
    """Calculates windspeed in m/s
    
    Parameters
    ----------
    obj : xarray dataset
        Model/obs data
    varmap : dictionary
        Dictionary defining the data column names to use
        Please provide variable names for "u_comp" in m/s, "v_comp" in m/s
    output_key : string
        String to name the new windspeed ouput
    Returns
    -------
    xarray dataset
        Xarray dataset with applied calculation as a new data array
        
    """
    
    # grab variable names from the yaml 
    u_key = varmap["u_comp"] 
    v_key = varmap["v_comp"] 

    u = obj[u_key]
    v = obj[v_key]

    wspd = (metpy.calc.wind_speed(
        u * units("m/s"),
        v * units("m/s")
        )).metpy.convert_units("m/s").metpy.dequantify()
    #Fancy dequantify method puts back in standard xarray format

    obj[output_key] = wspd
        
    return obj

# calc wind direction
def wdir(obj, varmap = None, output_key = "winddir"):
    """Calculates wind direction in degrees
    
    Parameters
    ----------
    obj : xarray dataset
        Model/obs data
    varmap : dictionary
        Dictionary defining the data column names to use
        Please provide variable names for "u_comp" in m/s, "v_comp" in m/s
    output_key : string
        String to name the new wind direction ouput
    Returns
    -------
    xarray dataset
        Xarray dataset with applied calculation as a new data array
        
    """
    
    # grab variable names from the yaml 
    u_key = varmap["u_comp"] 
    v_key = varmap["v_comp"] 
    
    u = obj[u_key].compute()
    v = obj[v_key].compute()
    #Unfortunately, wind_direction does not work with dask in metpy.
    #Maybe find another solution as we move to optimize memory usage, 
    #but for now just compute these.

    wdir = (metpy.calc.wind_direction(
        u * units("m/s"),
        v * units("m/s")
        )).metpy.convert_units("degree").metpy.dequantify()
    #Fancy dequantify method puts back in standard xarray format

    obj[output_key] = wdir
        
    return obj

# calc potential temperature
def ptemp(obj, varmap=None, output_key="ptemp"):
    """Calculates potential temperature in K
    
    Parameters
    ----------
    obj : xarray dataset
        Model/obs data
    varmap : dictionary
        Dictionary defining the data column names to use
        Please provide variable names for "pres_calc" in Pa, "temp_calc" in K
    output_key : string
        String to name the new potential temperature ouput
    Returns
    -------
    xarray dataset
        Xarray dataset with applied calculation as a new data array
        
    """

    pressure_key = varmap['pres_calc'] if varmap and 'pres_calc' in varmap else 'surfpres_pa'
    temperature_key = varmap['temp_calc'] if varmap and 'temp_calc' in varmap else 'temperature_k'
        
    pres = obj[pressure_key]
    temp = obj[temperature_key]
    
    ptmp = (metpy.calc.potential_temperature(
        pres * units.Pa,
        temp * units("K")
    )).metpy.convert_units("K").metpy.dequantify()
    #Fancy dequantify method puts back in standard xarray format

    obj[output_key] = ptmp
        
    return obj

