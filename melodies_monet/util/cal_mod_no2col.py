# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

# calculate WRF-Chem NO2 trop. columns, for further pair with satellite swath data
# contact: meng.li.atm@gmail.com
#

import xesmf as xe
import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime

def mod_to_overpasstime(modobj,opass_tms):
    '''
    Interpolate model to satellite overpass time.

    Parameters
    ----------
    modobj : xarray, model data
    opass_tms : datetime64, satellite overpass local time

    Output
    ------
    outmod : revised model data at local overpass time
    '''
    nst, = opass_tms.shape
    nmt, = modobj.time.shape
    ny,nx = modobj.longitude.shape
    
    # Determine local time offset
    local_utc_offset = np.zeros([ny,nx],dtype='timedelta64[ns]')
    # pandas timedelta calculation doesn't work on ndarrays
    for xi in np.arange(nx):
        local_utc_offset[:,xi] = pd.to_timedelta((modobj['longitude'].isel(x=xi)/15).astype(np.int64),unit='h')
    
    # initialize local time as variable
    modobj['localtime'] = (['time','y','x'],np.zeros([nt,ny,nx],dtype='datetime64[ns]'))
    # fill
    for ti in np.arange(nt):
        modobj['localtime'][ti] = modobj['time'][ti].data + local_utc_offset

    # initalize new model object with satellite datetimes
    outmod = []

    for ti in np.arange(nst):
        # Apply filter to select model data within +/- 1 output time step of the overpass time
        tempmod = modobj.where(np.abs(modobj['localtime'] - opass_tms[ti].to_datetime64()) < (modobj.time[1] - modobj.time[0]))
        
        # determine factors for linear interpolation in time
        tfac = 1 - (np.abs(tempmod['localtime'] - opass_tms[ti].to_datetime64())/(modobj.time[1] - modobj.time[0]))
        tempmod = tempmod.drop_vars('localtime')
        # Carry out time interpolation
        ## Note regarding current behavior: will only carry out time interpolation if at least 2 model timesteps
        outmod.append((tfac*tempmod).sum(dim='time', min_count=2,keep_attrs=True))
    outmod = xr.merge(outmod)
    outmod['time'] = (['time'],opass_tms)
    return outmod

def cal_model_no2columns(modobj):

    """
    Calcuate model (WRF-Chem) NO2 columns for each layer, to pair with satellite data
    
    Parameters
    ------
    modobj         : model data
  
    Output
    ------
    modobj        : revised model data with 'no2col' and 'localtime' added

    """

    # calculate the no2 tropospheric vertical columns and pressure from wrf-chem
    # update, mli, to be coordinated with monetio
    no2    = modobj['no2']
    tb     = modobj['temperature_k']
    pb2    = modobj['pres_pa_mid'] # time,z,y,x
    dzh    = modobj['dz_m'] # time,z,y,x, the model layer thickness 
    time   = modobj.coords['time']
    modlon = modobj.coords['longitude']


    # presure: base state + PB (KSMP)
    nt, nz, ny, nx = no2.shape
    #pb2  = np.zeros([nt, nz, ny, nx],dtype=np.float)
    #pb2  = pdata + pb

    # convert the perturbation potential temperature (from 300K reference) to temp
    #tb = np.zeros([nt, nz, ny, nx],dtype=np.float)
    #tb =(300.0+tdata)*((pb2/1.0e5)**0.286)


    # --- initialize arrays
    # no2 columns for each layer
    no2col     = np.zeros([nt, nz, ny, nx], dtype = np.float32)
    # temporary array
    value      = np.zeros([nt, ny, nx], dtype = np.float32)

    # average between 13:00 and 14:00 localtime
    localtm    = np.zeros([nt,ny,nx], dtype='datetime64[s]')
    tdlist     = np.zeros([ny], dtype=np.int32)
    tdlt       = np.zeros([ny, nx], dtype='timedelta64[ms]')

    for xx in range(nx):
         tdlist[:]  = np.array(modlon[:,xx]/15.0).astype(int)
         tdlt[:,xx] = pd.TimedeltaIndex(tdlist, 'h')

    for tt in range(nt):
        localtm[tt,:,:] = pd.Timestamp(time.values[tt]) + tdlt[:,:]

    # --- calculate NO2 columns by layers
    # convert to ppm
    no2 = no2 / 1000.0
    for vl in range(nz):
        ad = pb2[:,vl,:,:] * (28.97e-3)/(8.314*tb[:,vl,:,:])
        #zh = ((ph[:,vl+1,:,:] + phb[:,vl+1,:,:]) - (ph[:,vl,:,:]+phb[:,vl,:,:]))/9.81
        value[:,:,:]= no2[:,vl,:,:]*dzh[:,vl,:,:]*6.022e23/(28.97e-3)*1e-10*ad[:,:,:] # timex y x x
        no2col[:,vl,:,:] = value[:,:,:]

    # add to model
    #modobj['PB2'] = xr.DataArray(pb2,dims=["time", "z", "y","x"]) # change from "time" to "t"
    modobj['localtime'] = xr.DataArray(localtm, dims=["t","y", "x"])
    modobj['no2col'] = xr.DataArray(no2col,dims=["t", "z", "y","x"])

    return modobj
