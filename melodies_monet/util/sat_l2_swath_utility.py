# SPDX-License-Identifier: Apache-2.0
#

# read all swath data for the time range
# developed for TROPOMI Level2 NO2
#

import numpy as np
import xarray as xr
from datetime import datetime
import xesmf as xe

import logging
numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)

def trp_interp_swatogrd(obsobj, modobj,no2varname='no2'):

    """
    interpolate sat swath to model grid
    
    Parameters
    ------
    obsobj  : satellite swath data
    modobj  : model data (with no2 col calculated)
    
    Output
    ------
    no2_modgrid_avg: Regridded satellite data at model grids for all datetime

    """
    
    # model grids attributes
    nmodt, nz, ny, nx  = modobj[f'{no2varname}_col'].shape # time, z, y, x, no2 columns at molec cm^-2
    
    time   = [datetime.strptime(x,'%Y-%m-%d') for x in obsobj.keys()]
    nobstime  = len(list(obsobj.keys()))

    # daily averaged sat data at model grids
    no2_modgrid_avg=xr.Dataset(data_vars = {
        'nitrogendioxide_tropospheric_column':(["time", "x", "y"],
                                                np.full([nobstime, ny, nx], np.nan, dtype=np.float32)),
        f'{no2varname}trpcol':(["time", "x", "y"], np.full([nobstime, ny, nx], np.nan, dtype=np.float32))
            },
        coords = dict(
            time=time,
            longitude=(["x", "y"], modobj.coords['longitude'].values),
            latitude=(["x", "y"], modobj.coords['latitude'].values)),
        attrs=dict(description="daily tropomi data at model grids"),)

    for nd in range(nobstime):
        days = list(obsobj.keys())[nd]
        # --- model
        # get model no2 trop. columns at 13:00 - 14:00 localtime
        modobj_tm = modobj.sel(time=days.strfime('%Y-%d-%m'))
        
        # intermediate need: model NO2 partial columns for day
        # no2col_satm = np.nanmean(modobj_tm['no2col'].values, axis = 0)
        
        # sum up tropopause
        if 'pres_pa_trop' in list(modobj.keys()):
            no2_modgrid_avg[f'{no2varname}trpcol'][nd, :,:] = modobj_tm[f'{no2varname}_col'].where(modobj_tm['pres_pa_mid'] >= modobj_tm['pres_pa_trop']).sum(dim='z').values.squeeze()

        else:
            print('Caution: model tropospheric NO2 column was calculated assuming the model top is the tropopause')
            no2_modgrid_avg[f'{no2varname}trpcol'][nd, :,:] = modobj_tm[f'{no2varname}_col'].sum(dim='z').values.squeeze()
            
        # --- TROPOMI
        # number of swath
        nswath = len(obsobj[days])

        # intermediate array for all swaths
        no2_modgrid_all = np.zeros([ny, nx, nswath], dtype=np.float64)

        for ns in range(nswath):
            satlon = obsobj[days][ns]['lon']
            satlat = obsobj[days][ns]['lat']
            satno2 = obsobj[days][ns]['nitrogendioxide_tropospheric_column']

            # regridding from swath grid to model grids
            grid_in = {'lon':satlon.values, 'lat':satlat.values}

            regridder = xe.Regridder(grid_in, no2_modgrid_avg[['lat','lon']],'bilinear',ignore_degenerate=True,reuse_weights=False)
            
            # regridded no2 trop. columns
            no2_modgrid = regridder(satno2) # , keep_attrs=True
            print('Done with TROPOMI regridding', days, ns)

            #regridder.destroy()
            del regridder
 
            no2_modgrid_all[:,:,ns] = no2_modgrid
            print(' no2 satellite:', np.nanmin(no2_modgrid), np.nanmax(no2_modgrid))

        # daily averaged no2 trop. columns at model grids
        no2_modgrid_avg['nitrogendioxide_tropospheric_column'][nd,:,:] = np.nanmean(np.where(no2_modgrid_all > 0.0, no2_modgrid_all, np.nan), axis=2)

    del(modobj)
    del(obsobj)

    return no2_modgrid_avg


def trp_interp_swatogrd_ak(obsobj, modobj,no2varname='no2'):

    """
    interpolate sat swath to model grid applied with averaging kernel
    
    Parameters
    ------
    obsobj  : satellite swath data
    modobj  : model data (with no2 col calculated)
    
    Output
    ------
    no2_modgrid_avg: Regridded satellite data at model grids for all datetime

    """

    # model grids attributes
    nmodt, nz, ny, nx  = modobj[f'{no2varname}_col'].shape # time, z, y, x, no2 columns at molec cm^-2
    
    time   = [datetime.strptime(x,'%Y-%m-%d') for x in obsobj.keys()]
    nobstime  = len(list(obsobj.keys()))

    # daily averaged sat data at model grids
    no2_modgrid_avg=xr.Dataset(data_vars = {
        'nitrogendioxide_tropospheric_column':(["time", "x", "y"],
                                                np.full([nobstime, ny, nx], np.nan, dtype=np.float32)),
        f'{no2varname}trpcol':(["time", "x", "y"],np.full([nobstime, ny, nx], np.nan, dtype=np.float32))
            },
        coords = dict(
            time=time,
            longitude=(["x", "y"], modobj.coords['longitude'].values),
            latitude=(["x", "y"], modobj.coords['latitude'].values)),
        attrs=dict(description="daily tropomi data at model grids"),)

    # tmpvalue = np.zeros([ny, nx], dtype = np.float64)

    # loop over all days
    for nd in range(nobstime):

        days = time[nd].strftime('%Y-%m-%d')
        # --- model ---
        # get model no2 trop. columns at 13:00 - 14:00 localtime
        try:
            modobj_tm = modobj.sel(time=days)
        except KeyError:
            print(days)
            print('Satellite data was outside available model times')
            continue
        #modobj_tm = modobj.sel(time=days)
        # no2col_satm = modobj_tm[f'{no2varname}_col'].mean(dim='time')
              
        # sum up tropopause, needs to be revised to tropopause
        if 'pres_pa_trop' in list(modobj.keys()):
            no2_modgrid_avg[f'{no2varname}trpcol'][nd, :,:] = modobj_tm[f'{no2varname}_col'].where(modobj_tm['pres_pa_mid'] >= modobj_tm['pres_pa_trop']).sum(dim='z').values.squeeze()

        else:
            print('Caution: model tropospheric NO2 column was calculated assuming the model top is the tropopause')
            no2_modgrid_avg[f'{no2varname}trpcol'][nd, :,:] = modobj_tm[f'{no2varname}_col'].sum(dim='z').values.squeeze()
        # --- tropomi ---
        # number of swath
        nswath = len(obsobj[days])

        # array for all swaths
        no2_modgrid_all = np.zeros([ny, nx, nswath], dtype=np.float32)

        for ns in range(nswath):
            working_swath = obsobj[days][ns]     

            grid_sat = {'lon':working_swath['lon'].values, 'lat':working_swath['lat'].values}
            grid_mod= {'lon':modobj.coords['longitude'].values, 'lat':modobj.coords['latitude'].values}


            nysat, nxsat, nzsat = working_swath['averaging_kernel'].shape

            # regridding from model grid to sat grid
            regridder_ms = xe.Regridder(grid_mod, grid_sat,'bilinear',ignore_degenerate=True,reuse_weights=False)
            
            # force model data to put z dimension last for pressure and no2 partial columns
            mod_pres_no2 = modobj_tm[['pres_pa_mid',f'{no2varname}_col']].mean(dim='time')#.transpose('y','x','z')
            #print(mod_pres_no2['no2col'].shape)
            # regridding for model pressure, and no2 vertical columns
            mod_rgd_sat = regridder_ms(mod_pres_no2)
            mod_rgd_sat = mod_rgd_sat.transpose('y','x','z')
            # convert from aks to trop.aks
            working_swath['averaging_kernel'] = working_swath['averaging_kernel'] * working_swath['air_mass_factor_total'] / working_swath['air_mass_factor_troposphere']
            # calculate the revised tamf_mod, and ratio = tamf_mod / tamf_org
            ratio = cal_amf_wrfchem(working_swath['averaging_kernel'], mod_rgd_sat['pres_pa_mid'].values, working_swath['preslev'], working_swath['troppres'], mod_rgd_sat[f'{no2varname}_col'].values,
                                    working_swath['air_mass_factor_troposphere'], grid_sat['lon'], grid_sat['lat'], grid_mod['lon'], grid_mod['lat'])

            # averaing kernel applied done
            satno2 = working_swath['nitrogendioxide_tropospheric_column'] * ratio 

            # regridding from swath grid to model grids
            regridder = xe.Regridder(grid_sat, grid_mod,'bilinear',ignore_degenerate=True,reuse_weights=False)

            # regridded no2 trop. columns
            no2_modgrid = regridder(satno2, keep_attrs=True)
            no2_modgrid_all[:,:,ns] = no2_modgrid[:,:]

        # daily averaged no2 trop. columns at model grids
        no2_modgrid_avg['nitrogendioxide_tropospheric_column'][nd,:,:] = np.nanmean(np.where(no2_modgrid_all > 0.0, no2_modgrid_all, np.nan), axis=2)

    return no2_modgrid_avg


def cal_amf_wrfchem(scatw, wrfpreslayer, tpreslev, troppres, wrfno2layer_molec, tamf_org, satlon, satlat, modlon, modlat):
    from scipy import interpolate

    nsaty, nsatx, nz    = wrfpreslayer.shape
    nsatz, nsaty, nsatx = tpreslev.shape # mli, update to new dimension


    nume             = np.zeros([nsaty, nsatx], dtype=np.float32)
    deno             = np.zeros([nsaty, nsatx], dtype=np.float32)
    amf_wrfchem      = np.zeros([nsaty, nsatx], dtype=np.float32)
    amf_wrfchem[:,:] = np.nan
    wrfavk           = np.zeros([nsaty, nsatx, nz], dtype = np.float32)
    wrfavk[:,:,:]    = np.nan
    wrfavk_scl       = np.zeros([nsaty, nsatx], dtype=np.float32) 
    preminus         = np.zeros([nsaty, nsatx], dtype=np.float32)
    wrfpreslayer_slc = np.zeros([nsaty, nsatx], dtype=np.float32)
    tmpvalue_sat     = np.zeros([nsaty, nsatx], dtype=np.float32)
    tmpvalue_mod     = np.zeros([nsaty, nsatx], dtype=np.float32)
    
    
    # set the surface pressure to wrf one
    tpreslev[0,:,:] = wrfpreslayer[:,:,0] 

    # relationship between pressure to avk
    tpreslev = tpreslev.values 
    scatw    = scatw.values
    wrfpreslayer = np.where((wrfpreslayer <=0.0), np.nan, wrfpreslayer)

    # shrink the satellite domain to WRF
    lb = np.where( (satlon >= np.nanmin(modlon)) & (satlon <= np.nanmax(modlon)) 
        & (satlat >= np.nanmin(modlat)) & (satlat <= np.nanmax(modlat)))

    vertical_pres = []
    vertical_scatw = []
    vertical_wrfp = []
    
    if len(lb[0]) == 0:
        print('Caution: There are no observations within the model domain')
    for llb in range(len(lb[0])):
        yy = lb[0][llb]
        xx = lb[1][llb]
        vertical_pres = tpreslev[:,yy,xx] # mli, update to new dimension
        vertical_scatw = scatw[yy,xx,:]
        vertical_wrfp = wrfpreslayer[yy,xx,:]
        f = interpolate.interp1d(np.log10(vertical_pres[:]),vertical_scatw[:], fill_value="extrapolate")# relationship between pressure to avk
        wrfavk[yy,xx,:] = f(np.log10(vertical_wrfp[:])) #wrf-chem averaging kernel

    for l in range(nz-1):  # noqa: E741
        # check if it's within tropopause
        preminus[:,:]         = wrfpreslayer[:,:,l] - troppres[:,:]

        # wrfpressure and wrfavk
        wrfpreslayer_slc[:,:] = wrfpreslayer[:,:,l]
        wrfavk_scl[:,:]       = wrfavk[:,:,l]

        ind_ak = np.where(np.isinf(wrfavk_scl) | (wrfavk_scl <= 0.0))
        # use the upper level ak 
        if (ind_ak[0].size >= 1):
            tmpvalue_sat[:,:]  = wrfavk[:,:,l+1]
            wrfavk_scl[ind_ak] = tmpvalue_sat[ind_ak]

        ind = np.where(preminus >= 0.0)
        # within tropopause
        if (ind[0].size >= 1):
            # select grids that this level is within tropopause
            tmpvalue_mod[:,:]  = wrfno2layer_molec[:,:,l]
            nume[ind] += wrfavk_scl[ind]*tmpvalue_mod[ind]
            deno[ind] += tmpvalue_mod[ind]
        else:
            break
            
    # tropospheric amf calculated based on model profile and TROPOMI averaging kernel
    amf_wrfchem = nume / deno * tamf_org

    # ratio
    ratio = tamf_org / amf_wrfchem 

    # exclude nan
    ratio = np.where(np.isnan(ratio), 1.0, ratio)

    print('Done with Averaging Kernel revision,', 'factor min:',np.nanmin(ratio), 'max:',np.nanmax(ratio)) 

    return ratio 

# Conservative / unstructured TROPOMI NO2 operator (model -> TROPOMI column)
# Reuses the model-agnostic regrid() (regrid_util) + uxarray helpers

def _to_molmol(da):
    """Return model mixing ratio in mol/mol, deciding by VALUE MAGNITUDE.

    The CAM-unstructured reader scales NO2 to ppbV (~0.01-300) but can leave
    the units attribute as 'mol/mol' (stale-attr bug). mol/mol NO2 is
    ~1e-9-1e-7, so a max above ~1e-3 means the values are really ppbV.
    """
    mx = float(np.nanmax(np.abs(da.values)))
    return da * 1e-9 if mx > 1e-3 else da

def interp_vertical_mod2tropomi(obsobj, modobj_swath, variables=("NO2",)):
    """Interpolate model mixing ratio (intensive) onto TROPOMI's layers.
    
    interpolate the model *concentration* in log-pressure onto those layers; the partial
    column is integrated afterward (in :func:`apply_weights_mod2tropomi_no2`)
    using TROPOMI's own layer thickness 
    
    Parameters
    ----------
    obsobj : xr.Dataset
        One TROPOMI granule (time squeezed) with ``pres_pa_mid`` (z, y, x).
    modobj_swath : xr.Dataset
        Model already regridded to the swath pixels, with ``pres_pa_mid``
        (z_mod, y, x) and the requested ``variables``.

    Returns
    -------
    xr.Dataset
        ``variables`` interpolated onto TROPOMI's z layers, dims (z, y, x).
    """
    
    p_mod = np.asarray(modobj_swath["pres_pa_mid"].transpose("z", "y", "x").values)
    p_trop = np.asarray(obsobj["pres_pa_mid"].transpose("z", "y", "x").values)
    nz_t, ny, nx = p_trop.shape
    out = xr.Dataset()
    for var in list(variables):
        src = np.asarray(modobj_swath[var].transpose("z", "y", "x").values)
        dest = np.full((nz_t, ny, nx), np.nan, dtype=float)
        for j in range(ny):
            for i in range(nx):
                mp = p_mod[:, j, i]; mc = src[:, j, i]; tp = p_trop[:, j, i]
                good = np.isfinite(mp) & np.isfinite(mc)
                if good.sum() < 2 or not np.isfinite(tp).any():
                    continue
                order = np.argsort(mp[good])
                dest[:, j, i] = np.interp(
                    np.log10(tp), np.log10(mp[good][order]), mc[good][order],
                    left=np.nan, right=np.nan,
                )
        out[var] = xr.DataArray(dest, dims=("z", "y", "x"))
    return out

# regrid to lat lon 
# In the future, this might need a seperate util file that enables it to be generalizeable to other sat products 
_DEG_PER_M = 1.0 / 111320.0   # degrees latitude per metre (mean Earth)

def _model_lonlat_extent(modobj, pad=0.0):
    """(lonmin, lonmax, latmin, latmax) of the model domain, optionally padded"""
    mlon = np.asarray(modobj["longitude"].values)
    mlat = np.asarray(modobj["latitude"].values)
    return (float(np.nanmin(mlon)) - pad, float(np.nanmax(mlon)) + pad,
            float(np.nanmin(mlat)) - pad, float(np.nanmax(mlat)) + pad)


def _swath2latlon(swath, data_vars, res, extent, units = "deg", method="radius_mean"):
    """Regrid swath-paired fields (y, x) onto a regular lat/lon grid 

    use radius mean averaging 
    
    returns a Dataset on dims (lat, lon).
    """
    from melodies_monet.util.regrid_util import regrid

    lonmin, lonmax, latmin, latmax = extent
    if units in ("km", "m"):
        res_m = float(res) * 1000.0 if units == "km" else float(res)
        dlat = res_m * _DEG_PER_M
        midlat = 0.5 * (latmin + latmax)
        dlon = res_m * _DEG_PER_M / max(np.cos(np.deg2rad(midlat)), 0.1)
    else:  # degrees
        dlat = dlon = float(res)

    tlon1 = np.arange(lonmin, lonmax + dlon, dlon)
    tlat1 = np.arange(latmin, latmax + dlat, dlat)
    tlon2, tlat2 = np.meshgrid(tlon1, tlat1)        # (lat, lon)

    # handle the lat lon naming between satellites 
    # tropomi carry long / lat on (y,x) where tempo does lon / lat on (x,y)
    # rename to longintude latitude 
    ren = {}
    if "longitude" not in swath.variables and "lon" in swath.variables:
        ren["lon"] = "longitude"
    if "latitude" not in swath.variables and "lat" in swath.variables:
        ren["lat"] = "latitude"
    if ren:
        swath = swath.rename(ren)
    _to_coord = [c for c in ("longitude", "latitude") if c in swath.data_vars]
    if _to_coord:
        swath = swath.set_coords(_to_coord)

    # insert a conservative regriding option 
    # Builds the swath source mesh from corner bounds and a rectilinear target
    # mesh, then mesh-to-mesh conservative regrid
    # Fall back to radius_mean if bounds are missing
    
    if method in _CONSERVATIVE and "longitude_bounds" in swath.variables:
        try:
            import uxarray as ux

            swath_grid, _ = _tropomi_swath_mesh(swath)
            _hd = tuple(swath["longitude"].dims)
            _flat = (
                swath[data_vars].stack(n_face=_hd).reset_index("n_face", drop=True)
            )
            _flat = _flat.drop_vars(
                [c for c in _flat.coords if "n_face" in _flat[c].dims], errors="ignore")
            _src = ux.UxDataset(_flat, uxgrid=swath_grid)
            _tgt = ux.Grid.from_structured(lon=tlon1, lat=tlat1)
            _out = regrid(_src, method=method, target_grid=_tgt)
            _out = _out.where(_out != 0)        # empty cells go to  NaN, not 0

            _nlat, _nlon = tlat1.size, tlon1.size
            _ncell = _nlat * _nlon
            _fd = next((d for d in _out.dims if _out.sizes.get(d) == _ncell), None)
            if _fd is None:
                raise ValueError("conservative target face dim not found after regrid")
            _res = xr.Dataset()
            for v in _out.data_vars:
                da = _out[v]
                if _fd not in da.dims:
                    continue
    
    hdims = list(swath["longitude"].dims)            # (y, x) or (x, y)
    
    flat = (
        swath[data_vars]
        .stack(pixel=hdims)
        .reset_index("pixel", drop=True)
    )
    
    # search radius ~ cell size, but at least ~5 km (sensor footprint) so a
    # finer-than-sensor grid fills from the nearest pixel instead of going empty.
    radius = max(float(dlat), 0.05)
    out = regrid(flat, target={"lon": tlon2, "lat": tlat2},
                 method="radius_mean", radius=radius, target_dims=("lat", "lon"))
    out = out.where(out != 0)        # empty cells -> NaN, not 0

    out = out.assign_coords(lon=("lon", tlon1), lat=("lat", tlat1))

    # want to make sure these regrided lat lon pairs can just run through the existing plotting 
    return out.rename({"lat": "y", "lon": "x"})

def _swath_to_target(swath, modobj, method, data_vars, target, res, extent, units="deg" ):
    """Regrid the paired swath onto the requested target space

    if target = model, use model's native unstructured grid via tropomi_swath2mod

    if target = obs, use a regular lat lon grid via _swath2latlon 

    """
    if target == "model":
        return _tropomi_swath2mod(swath, modobj, method, data_vars)
    if target == "obs":
        return _swath2latlon(swath, data_vars, res, extent, units=units, method=method)
    raise ValueError(f"regrid_target {target!r} not understood; use 'model' or 'obs'.")

def apply_weights_mod2tropomi_no2(obsobj, modobj_on_tropomi_layers, species="NO2"):
    """Apply the TROPOMI averaging kernel to a model NO2 profile.

    Mirrors :func:`apply_weights_mod2tempo_no2` but uses TROPOMI's
    averaging kernel instead of scattering weights:

      AK_trop = averaging_kernel * (amf_total / amf_troposphere)
      subcol  = vmr * (dp/g) * (NA/M_air) / 1e4         # molec/cm2 per layer
      VCD     = sum over tropospheric layers (p >= tropopause) of AK_trop*subcol

    Parameters
    ----------
    obsobj : xr.Dataset
        One TROPOMI granule (time squeezed): averaging_kernel (z,y,x),
        air_mass_factor_troposphere/_total (y,x), pres_pa_mid (z,y,x),
        pres_pa_int (z_stagg,y,x), tm5_tropopause_pressure (y,x).
    modobj_on_tropomi_layers : xr.Dataset
        Model ``species`` mixing ratio on TROPOMI's z layers (z,y,x), from
        :func:`interp_vertical_mod2tropomi`. Units auto-detected (ppbV/mol/mol).

    Returns
    -------
    xr.DataArray
        Model NO2 tropospheric column with the AK applied, molec/cm2, (y, x).
    """
    g, M_air, NA = 9.80665, 0.0289644, 6.022e23

    vmr = _to_molmol(modobj_on_tropomi_layers[species]).transpose("z", "y", "x")

    pint = obsobj["pres_pa_int"].transpose("z_stagg", "y", "x")
    dp = np.abs(
        pint.isel(z_stagg=slice(0, -1)).values
        - pint.isel(z_stagg=slice(1, None)).values
    )
    dp = xr.DataArray(dp, dims=("z", "y", "x"))
    subcol = vmr * dp * (NA / (g * M_air) / 1e4)            # molec/cm2 per layer

    ak = obsobj["averaging_kernel"].transpose("z", "y", "x")
    ak_trop = ak * (obsobj["air_mass_factor_total"]
                    / obsobj["air_mass_factor_troposphere"])
    trop = (obsobj["pres_pa_mid"].transpose("z", "y", "x")
            >= obsobj["tm5_tropopause_pressure"])

    vcd = (ak_trop * subcol).where(trop).sum("z", skipna=True)
    vcd = vcd.where(np.isfinite(vmr.isel(z=0)))
    vcd.attrs = {
        "units": "molecules/cm2",
        "description": "model NO2 tropospheric column after applying TROPOMI averaging kernel",
        "history": "Created by MELODIES-MONET, apply_weights_mod2tropomi_no2",
    }
    return vcd.where(np.isfinite(vcd))

# Orchestrator: conservative/unstructured TROPOMI NO2 pairing.
# Mirrors the TEMPO regrid_and_apply_weights flow but with the TROPOMI
# averaging-kernel operator. Reuses only the shared, product-agnostic regrid
# primitives (regrid_util.regrid, uxarray_util.open_uxgrid /
# uxgrid_from_corner_bounds) -- no dependency on the TEMPO utility.

_TROPOMI_NO2_VAR = "nitrogendioxide_tropospheric_column"
_MOL_M2_TO_MOLEC_CM2 = 6.02214e19
_CONSERVATIVE = ("conservative", "conservative_normed")

def _tropomi_swath_mesh(o):
    """Build a uxarray Grid from a TROPOMI granule's pixel corner bounds.

    Returns (grid, (ny, nx)). Flatten order is row-major over (y, x), matching
    how the n_face result is reshaped back.
    """
    from melodies_monet.util.uxarray_util import uxgrid_from_corner_bounds

    olon = np.asarray(o["longitude"].values)  # (y, x)
    ny, nx = olon.shape
    clon = np.asarray(o["longitude_bounds"].values).reshape(ny * nx, -1)
    clat = np.asarray(o["latitude_bounds"].values).reshape(ny * nx, -1)
    return uxgrid_from_corner_bounds(clon, clat), (ny, nx)


def _mod2tropomi_swath(modobj, o, method, mod_vars, grid_file):
    """Regrid model fields onto the TROPOMI swath pixels (y, x).

    conservative -> mesh-to-mesh via xregrid (model SCRIP/MPAS mesh -> swath
    cells from bounds); nearest/radius -> cKDTree at pixel centers. Returns a
    Dataset with each requested var on (..., y, x).
    """
    from melodies_monet.util.regrid_util import regrid

    msrc = modobj[[v for v in mod_vars if v in modobj.variables]]
    olon = np.asarray(o["longitude"].values)
    olat = np.asarray(o["latitude"].values)
    ny, nx = olon.shape

    if method in _CONSERVATIVE:
        swath_grid, _ = _tropomi_swath_mesh(o)
        out = regrid(msrc, method=method, src_grid=grid_file, target_grid=swath_grid)
        nface = ny * nx
        res = xr.Dataset()
        for v in out.data_vars:
            da = out[v]
            fd = next((d for d in da.dims if da.sizes[d] == nface), None)
            if fd is None:
                continue
            t = da.transpose(..., fd)
            arr = np.asarray(t.values).reshape(t.shape[:-1] + (ny, nx))
            res[v] = xr.DataArray(arr, dims=t.dims[:-1] + ("y", "x"),
                                  attrs=dict(da.attrs))
        res = res.assign_coords(longitude=(("y", "x"), olon),
                                latitude=(("y", "x"), olat))
        return res

    # nearest family
    out = regrid(msrc, target={"lon": olon, "lat": olat},
                 method=method, target_dims=("y", "x"))
    return out


def _tropomi_swath2mod(swath, modobj, method, data_vars):
    """Regrid swath-paired fields (y, x) onto the unstructured model columns.

    conservative -> mesh-to-mesh via xregrid (swath cells -> model mesh);
    else -> cKDTree within-radius mean. Returns a Dataset on the model
    column dim with longitude/latitude attached.
    """
    import uxarray as ux
    from melodies_monet.util.regrid_util import regrid
    from melodies_monet.util.uxarray_util import open_uxgrid

    col_dim = modobj["longitude"].dims[0]
    mlon = np.asarray(modobj["longitude"].values).ravel()
    mlat = np.asarray(modobj["latitude"].values).ravel()
    grid_file = (modobj.attrs.get("mio_scrip_file")
                 or modobj.attrs.get("mio_mesh_file"))

    if method in _CONSERVATIVE and "longitude_bounds" in swath.variables:
        swath_grid, (ny, nx) = _tropomi_swath_mesh(swath)
        flat = (
            swath[data_vars]
            .stack(n_face=("y", "x"))
            .reset_index("n_face", drop=True)
        )
        flat = flat.drop_vars(
            [c for c in flat.coords if "n_face" in flat[c].dims], errors="ignore"
        )
        src = ux.UxDataset(flat, uxgrid=swath_grid)
        model_grid = open_uxgrid(grid_file)
        out = regrid(src, method=method, target_grid=model_grid)

        # Conservative regrid fills model cells with no swath overlap with exactly 0.
        out = out.where(out != 0)
        
        n_col = int(model_grid.n_face)
        d = next((dd for dd in out.dims if out.sizes[dd] == n_col), None)
        if d is not None and d != col_dim:
            out = out.rename({d: col_dim})
        return out.assign_coords({"longitude": (col_dim, mlon),
                                  "latitude": (col_dim, mlat)})

    # nearest / radius_mean fallback
    flat = (
        swath[data_vars]
        .stack(pixel=("y", "x"))
        .reset_index("pixel", drop=True)
        .set_coords(["longitude", "latitude"])
    )
    out = regrid(flat, target={"lon": mlon, "lat": mlat},
                 method="radius_mean", radius=0.1, target_dims=(col_dim,))

    # make sure to fill 0s with nans
    return out.where(out != 0)


def regrid_and_apply_weights_tropomi(obsobj, modobj, species=["NO2"],
                                     method="conservative", qa_min=0.75,
                                     regrid_target="model", obs_grid_res=0.1, obs_grid_units="deg", obs_grid_extent=None):
    
    """Pair an unstructured model with TROPOMI L2 NO2 (AK applied).

    For each granule: forward-regrid model to swat (obs), interpolate the model
    profile onto TROPOMI's layers, apply the averaging kernel, then
    regrid the AK'd model column AND the obs column back onto the model
    grid. Granules are concatenated along ``time``.

    Parameters
    ----------
    obsobj : dict[str, list[xr.Dataset]]
        Output of the generic TROPOMI reader (tropomi_l2.open_datasets):
        keyed by date, each value a list of orbit granules.
    modobj : xr.Dataset
        Unstructured model with longitude/latitude on its column dim,
        NO2 (ppbV or mol/mol), pres_pa_mid, and mio_scrip_file/mio_mesh_file.
    species : list[str]
        Model species name(s); species[0] is paired.
    method : str
        Regrid method (conservative recommended; nearest_s2d/radius_mean ok).

    Returns
    -------
    xr.Dataset
        Paired model + obs NO2 tropospheric columns (molec/cm2) on the model
        column dim, stacked along ``time`` (one entry per granule).
    """
    sp = species[0]
    grid_file = (modobj.attrs.get("mio_scrip_file")
                 or modobj.attrs.get("mio_mesh_file"))

    targets = [regrid_target] if isinstance(regrid_target, str) else list(regrid_target)
    extent = (tuple(obs_grid_extent) if obs_grid_extent
                  else _model_lonlat_extent(modobj))

    # Flatten dict[date -> list[granule]] to a flat granule list.
    granules = []
    for v in obsobj.values():
        granules.extend(v if isinstance(v, list) else [v])

    out_by = {t: [] for t in targets}
    for o in granules:
        if "time" in o.dims:
            o = o.squeeze("time", drop=False)
        
        # Granule overpass time for model matching. prefer "time_granule", which holds
        # the real per-measurement times, and use its mean. Fall back to the
        # "time" coord (midnight, start of day reference) only if "time_granule" is absent.
        
        if "time_granule" in o.variables:
            tg = np.asarray(o["time_granule"].values).ravel().astype("datetime64[ns]")
            tg = tg[~np.isnat(tg)]
            gtime = (
                np.array(tg.astype("int64").mean(), dtype="int64").astype("datetime64[ns]")
                if tg.size else None
            )
        else:
            gtime = o["time"].values if "time" in o.coords else None

        # Select the model to this granule's overpass time
        if "time" in modobj.dims and gtime is not None:
            tsel = gtime if np.ndim(gtime) == 0 else np.asarray(gtime).ravel()[0]

            mtimes = modobj["time"].values
            tmin, tmax = mtimes.min(), mtimes.max()
            if tsel < tmin:
                tsel = tmin
            elif tsel > tmax:
                tsel = tmax
                
            try:
                mod_t = modobj.interp(time=tsel)
            except Exception:
                mod_t = modobj.sel(time=tsel, method="nearest")
        elif "time" in modobj.dims:
            mod_t = modobj.isel(time=0)
        else:
            mod_t = modobj
            
        # forward model to swath (NO2 + model pressure for the vertical interp)
        on_swath = _mod2tropomi_swath(
            mod_t, o, method, ["NO2", "pres_pa_mid"], grid_file
        )
        # model vmr onto TROPOMI's 34 layers
        no2_t = interp_vertical_mod2tropomi(o, on_swath, ["NO2"])
        # averaging-kernel applied model column (molec/cm2), (y, x)
        model_col = apply_weights_mod2tropomi_no2(o, no2_t, "NO2")
        obs_col = o[_TROPOMI_NO2_VAR] * _MOL_M2_TO_MOLEC_CM2  # molec/cm2

        # qa val for no2
        if qa_min and "qa_value" in o.variables:
            qa = o["qa_value"]
            if float(np.nanmax(np.asarray(qa.values))) > 1.5:
                qa = qa / 100.0
            good = qa >= qa_min
            obs_col = obs_col.where(good)
            model_col = model_col.where(good)
                    
        swath = xr.Dataset(
            {sp: model_col, _TROPOMI_NO2_VAR: obs_col,
             "latitude_bounds": o["latitude_bounds"],
             "longitude_bounds": o["longitude_bounds"]},
            coords={"longitude": o["longitude"], "latitude": o["latitude"]},
        )

        for t in targets:
            on = _swath_to_target(swath, modobj, method, [sp, _TROPOMI_NO2_VAR],
                                  t, obs_grid_res, extent, units=obs_grid_units)
            if gtime is not None:
                tval = gtime if np.ndim(gtime) == 0 else np.asarray(gtime).ravel()[0]
                on = on.expand_dims(time=[np.datetime64(tval)])
            out_by[t].append(on)

    return {t: (xr.concat(lst, dim="time") if lst else xr.Dataset())
            for t, lst in out_by.items()}

# tropomi HCHO 
# single tropo AMF and total columng averaging kernel 
# use generic tropomi_l2 reader + shared regrid / interp 

def apply_weights_mod2tropomi_hcho(obsobj, modobj_on_tropomi_layers, species="CH2O"):
    """Apply the TROPOMI HCHO averaging kernel to a model formaldehyde profile.

    Returns the model HCHO column with the AK applied (molec/cm2), dims (y, x).
    """
    g, M_air, NA = 9.80665, 0.0289644, 6.022e23

    vmr = _to_molmol(modobj_on_tropomi_layers[species]).transpose("z", "y", "x")

    pint = obsobj["pres_pa_int"].transpose("z_stagg", "y", "x")
    dp = np.abs(
        pint.isel(z_stagg=slice(0, -1)).values
        - pint.isel(z_stagg=slice(1, None)).values
    )
    dp = xr.DataArray(dp, dims=("z", "y", "x"))
    subcol = vmr * dp * (NA / (g * M_air) / 1e4)            # molec/cm2 per layer

    # AK * AMF 
    ak = obsobj["averaging_kernel"].transpose("z", "y", "x")
    ak = ak * obsobj["formaldehyde_tropospheric_air_mass_factor"] # hopefully this isnt hardcoded otherwise will need to pull this in from YAML

    vcd = (ak * subcol).sum("z", skipna=True)
    vcd = vcd.where(np.isfinite(vmr.isel(z=0)))
    
    vcd.attrs = {
        "units": "molecules/cm2",
        "description": "model HCHO column after applying TROPOMI averaging kernel",
        "history": "Created by MELODIES-MONET, apply_weights_mod2tropomi_hcho",
    }
    return vcd.where(np.isfinite(vcd))

def regrid_and_apply_weights_tropomi_hcho(obsobj, modobj, species=["CH2O"],
                                          method="conservative", qa_min=0.5,
                                          regrid_target="model", obs_grid_res=0.1, obs_grid_units="deg", obs_grid_extent=None): 
    """
    Pair an unstructured model with TROPOMI L2 HCHO (AK applied).

    defaults to regridding to model space via "regrid_target"

    regrid_target="model" ; model space 
    regrid_target="obs" ; obs space
    
    """
    sp = species[0]
    grid_file = (modobj.attrs.get("mio_scrip_file")
                 or modobj.attrs.get("mio_mesh_file"))

    targets = [regrid_target] if isinstance(regrid_target, str) else list(regrid_target)
    extent = (tuple(obs_grid_extent) if obs_grid_extent
              else _model_lonlat_extent(modobj))

    granules = []
    for v in obsobj.values():
        granules.extend(v if isinstance(v, list) else [v])

    out_by = {t: [] for t in targets}
    for o in granules:
        if "time" in o.dims:
            o = o.squeeze("time", drop=False)

        if "time_granule" in o.variables:
            tg = np.asarray(o["time_granule"].values).ravel().astype("datetime64[ns]")
            tg = tg[~np.isnat(tg)]
            gtime = (
                np.array(tg.astype("int64").mean(), dtype="int64").astype("datetime64[ns]")
                if tg.size else None
            )
        else:
            gtime = o["time"].values if "time" in o.coords else None

        if "time" in modobj.dims and gtime is not None:
            tsel = gtime if np.ndim(gtime) == 0 else np.asarray(gtime).ravel()[0]
            mtimes = modobj["time"].values
            tmin, tmax = mtimes.min(), mtimes.max()
            if tsel < tmin:
                tsel = tmin
            elif tsel > tmax:
                tsel = tmax
            try:
                mod_t = modobj.interp(time=tsel)
            except Exception:
                mod_t = modobj.sel(time=tsel, method="nearest")
        elif "time" in modobj.dims:
            mod_t = modobj.isel(time=0)
        else:
            mod_t = modobj

        on_swath = _mod2tropomi_swath(mod_t, o, method, [sp, "pres_pa_mid"], grid_file)
        prof_t = interp_vertical_mod2tropomi(o, on_swath, [sp])
        model_col = apply_weights_mod2tropomi_hcho(o, prof_t, sp)
        obs_col = o["formaldehyde_tropospheric_vertical_column"] * _MOL_M2_TO_MOLEC_CM2     # molec/cm2

        # QA filter
        if qa_min and "qa_value" in o.variables:
            qa = o["qa_value"]
            if float(np.nanmax(np.asarray(qa.values))) > 1.5:
                qa = qa / 100.0
            good = qa >= qa_min
            obs_col = obs_col.where(good)
            model_col = model_col.where(good)
                    
        swath = xr.Dataset(
            {sp: model_col, "formaldehyde_tropospheric_vertical_column": obs_col,
             "latitude_bounds": o["latitude_bounds"],
             "longitude_bounds": o["longitude_bounds"]},
            coords={"longitude": o["longitude"], "latitude": o["latitude"]},
        )

        for t in targets:
            on = _swath_to_target(
                swath, modobj, method,
                [sp, "formaldehyde_tropospheric_vertical_column"],
                t, obs_grid_res, extent, units=obs_grid_units
            )
            if gtime is not None:
                tval = gtime if np.ndim(gtime) == 0 else np.asarray(gtime).ravel()[0]
                on = on.expand_dims(time=[np.datetime64(tval)])
            out_by[t].append(on)

    return {t: (xr.concat(lst, dim="time") if lst else xr.Dataset())
            for t, lst in out_by.items()}


# new satelite variables ncdump -h "$F" | grep -iE "group: (PRODUCT|DETAILED_RESULTS|INPUT_DATA)|carbonmonoxide_total_column|column_averaging_kernel|pressure_levels|qa_value|layer =|:units|:multiplication_factor"

# co 

_TROPOMI_CO_VAR = "carbonmonoxide_total_column"

def apply_weights_mod2tropomi_co(obsobj, modobj_on_tropomi_layers, species="CO"):
    """Apply the TROPOMI CO column averaging kernel to a model CO profile.

    """
    g, M_air, NA = 9.80665, 0.0289644, 6.022e23

    vmr = _to_molmol(modobj_on_tropomi_layers[species]).transpose("z", "y", "x")

    pint = obsobj["pres_pa_int"].transpose("z_stagg", "y", "x")
    dp = np.abs(
        pint.isel(z_stagg=slice(0, -1)).values
        - pint.isel(z_stagg=slice(1, None)).values
    )
    dp = xr.DataArray(dp, dims=("z", "y", "x"))
    subcol = vmr * dp * (NA / (g * M_air) / 1e4)            # molec/cm2 per layer

    ak = obsobj["column_averaging_kernel"].transpose("z", "y", "x")  # dimensionless 
    vcd = (ak * subcol).sum("z", skipna=True)
    vcd = vcd.where(np.isfinite(vmr.isel(z=0)))
    
    vcd.attrs = {
        "units": "molecules/cm2",
        "description": "model CO column after applying TROPOMI column averaging kernel",
        "history": "Created by MELODIES-MONET, apply_weights_mod2tropomi_co",
    }
    return vcd.where(np.isfinite(vcd))

def regrid_and_apply_weights_tropomi_co(obsobj, modobj, species=["CO"],
                                        method="conservative", qa_min=0.5,
                                        regrid_target="model", obs_grid_res=0.1, obs_grid_units="deg", obs_grid_extent=None):
    
    """Pair an unstructured model with TROPOMI L2 CO (column AK applied).

    Mirrors regrid_and_apply_weights_tropomi_hcho 

    defaults to regridding to model space via "regrid_target"

    regrid_target="model" ; model space 
    regrid_target="obs" ; obs space

    """
    sp = species[0]
    grid_file = (modobj.attrs.get("mio_scrip_file")
                 or modobj.attrs.get("mio_mesh_file"))

    targets = [regrid_target] if isinstance(regrid_target, str) else list(regrid_target)
    extent = (tuple(obs_grid_extent) if obs_grid_extent
              else _model_lonlat_extent(modobj))

    granules = []
    for v in obsobj.values():
        granules.extend(v if isinstance(v, list) else [v])

    out_by = {t: [] for t in targets}
    for o in granules:
        if "time" in o.dims:
            o = o.squeeze("time", drop=False)

        if "time_granule" in o.variables:
            tg = np.asarray(o["time_granule"].values).ravel().astype("datetime64[ns]")
            tg = tg[~np.isnat(tg)]
            gtime = (
                np.array(tg.astype("int64").mean(), dtype="int64").astype("datetime64[ns]")
                if tg.size else None
            )
        else:
            gtime = o["time"].values if "time" in o.coords else None

        if "time" in modobj.dims and gtime is not None:
            tsel = gtime if np.ndim(gtime) == 0 else np.asarray(gtime).ravel()[0]
            mtimes = modobj["time"].values
            tmin, tmax = mtimes.min(), mtimes.max()
            if tsel < tmin:
                tsel = tmin
            elif tsel > tmax:
                tsel = tmax
            try:
                mod_t = modobj.interp(time=tsel)
            except Exception:
                mod_t = modobj.sel(time=tsel, method="nearest")
        elif "time" in modobj.dims:
            mod_t = modobj.isel(time=0)
        else:
            mod_t = modobj

        on_swath = _mod2tropomi_swath(mod_t, o, method, [sp, "pres_pa_mid"], grid_file)
        prof_t = interp_vertical_mod2tropomi(o, on_swath, [sp])
        model_col = apply_weights_mod2tropomi_co(o, prof_t, sp)
        obs_col = o[_TROPOMI_CO_VAR] * _MOL_M2_TO_MOLEC_CM2     # molec/cm2

        # QA filter 
        if qa_min and "qa_value" in o.variables:
            qa = o["qa_value"]
            if float(np.nanmax(np.asarray(qa.values))) > 1.5:
                qa = qa / 100.0
            good = qa >= qa_min
            obs_col = obs_col.where(good)
            model_col = model_col.where(good)

        swath = xr.Dataset(
            {sp: model_col, _TROPOMI_CO_VAR: obs_col,
             "latitude_bounds": o["latitude_bounds"],
             "longitude_bounds": o["longitude_bounds"]},
            coords={"longitude": o["longitude"], "latitude": o["latitude"]},
        )

        for t in targets:
            on = _swath_to_target(swath, modobj, method, [sp, _TROPOMI_CO_VAR],
                                  t, obs_grid_res, extent, units=obs_grid_units)
            if gtime is not None:
                tval = gtime if np.ndim(gtime) == 0 else np.asarray(gtime).ravel()[0]
                on = on.expand_dims(time=[np.datetime64(tval)])
            out_by[t].append(on)

    return {t: (xr.concat(lst, dim="time") if lst else xr.Dataset())
            for t, lst in out_by.items()}

