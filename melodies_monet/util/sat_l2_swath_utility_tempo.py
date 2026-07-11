# SPDX-License-Identifier: Apache-2.0
#

# read all swath data for the time range
# developed for TEMPO Level2 NO2
#

"""Python utility for TEMPO use."""

import collections
import glob
import gc
import logging
import warnings

import numba
import numpy as np
import xarray as xr
import xesmf as xe

numba_logger = logging.getLogger("numba")
numba_logger.setLevel(logging.WARNING)

# Regridding methods valid for unstructured (cell-data) models such as
# CESM-SE / MUSICA. ESMF 'bilinear' and 'patch' interpolate from mesh
# NODES, but model output lives on cell centers (faces) -- DO NOT USE. 
# Conservative (face-to-face, mass-conserving) is recommended;
# the nearest family is the fast, non-conserving alternative.

_UNSTRUCT_CONSERVATIVE = ("conservative", "conservative_normed")
_UNSTRUCT_NEAREST = ("nearest_s2d", "nearest_d2s", "radius_mean")
_UNSTRUCT_SUPPORTED = _UNSTRUCT_CONSERVATIVE + _UNSTRUCT_NEAREST

def _unsupported_method_error(method, where):
    """ValueError explaining which regrid methods unstructured
    models support, and why bilinear/patch don't."""
    return ValueError(
        f"{where}: regrid_method={method!r} is not supported for unstructured "
        "(cell-data) models such as CESM-SE / MUSICA.\n"
        "ESMF 'bilinear' and 'patch' interpolate from mesh NODES, but the "
        "model output is on cell centers (faces), so those methods cannot be "
        "applied.\n"
        "Supported methods (set 'regrid_method:' in the obs YAML block):\n"
        "  - 'conservative'  (RECOMMENDED) true area-weighted, mass-conserving\n"
        "  - 'nearest_s2d' / 'nearest_d2s'  fast nearest-neighbor (not conserving)\n"
        "  - 'radius_mean'  within-radius average (swath->model direction)"
    )

def calc_grid_corners(ds, lat="latitude", lon="longitude"):
    """Adds latitude and longitude bounds inplace.
    If the grid is rectilinear, it should be quite precise.
    If it is curvilinear, is a rough estimate.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to which the latitude and longitude will be added.
    lat : str
        name of the lat variable.
    lon : str
        name of the lon variable.

    Returns
    -------
    None
    """
    try:
        import cf_xarray as cfxr
    except ImportError:
        raise ImportError("Calculating gridcell bounds requires cf_xarray. Please install")
    corners = ds[[lat, lon]].cf.add_bounds([lat, lon])
    ds["lat_b"] = cfxr.bounds_to_vertices(corners[f"{lat}_bounds"], "bounds", order=None)
    ds["lon_b"] = cfxr.bounds_to_vertices(corners[f"{lon}_bounds"], "bounds", order=None)
    return


def speedup_regridding(dset, variables="all"):
    """Makes modobj latitude and longitude C_contiguous, which speeds up regridding.
    It makes the changes inplace

    Parameters
    ----------
    dset : xr.Dataset
        Dataset containing latitude and longitude

    Returns
    -------
    None
    """
    if variables == "all":
        variables = dset.variables
    for v in variables:
        if not dset[v].values.flags["C_CONTIGUOUS"]:
            dtype = dset[v].dtype
            dset[v] = dset[v].astype(dtype, order="C")

def _nearest_mod2swath(modobj, obsobj):
    """Sample an unstructured model at swath pixel locations.

    Thin wrapper around :func:`melodies_monet.util.regrid.regrid` for the
    forward (model -> swath) direction: the 2-D swath ``(x, y)`` lon/lat
    are the targets, and the model column centers are the cKDTree source.
    Downstream swath code (vertical interp, apply weights) expects
    ``lon``/``lat`` (short names) on ``(x, y)``, so we rename the
    long-name coords ``regrid`` attaches back to the swath convention.
    
    """
    from melodies_monet.util.regrid_util import regrid

    out = regrid(
        modobj,
        target={
            "lon": np.asarray(obsobj["lon"].values),
            "lat": np.asarray(obsobj["lat"].values),
        },
        method="nearest_s2d",
        target_dims=("x", "y"),
    )

    return out.rename({"longitude": "lon", "latitude": "lat"})

def _swath_corner_bounds(obsobj):
    """Return (lon_bounds, lat_bounds) DataArrays for a TEMPO swath.

    Handles the naming variants: TEMPO V03 ``latitude_bounds`` /
    ``longitude_bounds``, generic ``lat_bounds`` / ``lon_bounds``, or the
    ``calc_grid_corners`` output ``lat_b`` / ``lon_b``.
    """
    for lon_n, lat_n in (
        ("longitude_bounds", "latitude_bounds"),
        ("lon_bounds", "lat_bounds"),
        ("lon_b", "lat_b"),
    ):
        if lon_n in obsobj.variables and lat_n in obsobj.variables:
            return obsobj[lon_n], obsobj[lat_n]
    return None, None

def _carry_swath_bounds(output, obs_src):
    """Merge swath corner bounds from ``obs_src`` into ``output`` if present.

    Lets the backward (swath -> model) conservative path rebuild the swath
    mesh later: the paired output otherwise drops the bounds.
    """
    lon_b, lat_b = _swath_corner_bounds(obs_src)
    if lon_b is None:
        return output
    extra = xr.Dataset({lon_b.name: lon_b, lat_b.name: lat_b})
    return xr.merge([output, extra], compat="override")


def _conservative_mod2swath(modobj, obsobj, method="conservative"):
    """Conservatively regrid an unstructured model onto a TEMPO swath.

    Builds a uxarray Grid from the swath pixel corner bounds
    (``longitude_bounds``/``latitude_bounds``, shape ``(x, y, corner)``),
    then uses xregrid mesh-to-mesh conservative regridding from the model's
    (depadded SCRIP) mesh onto the swath cells. Result is reshaped back to
    the swath ``(..., x, y)`` layout with ``lon``/``lat`` coords, matching
    the convention downstream vertical-interp / apply-weights code expects.
    """
    from melodies_monet.util.regrid_util import regrid
    from melodies_monet.util.uxarray_util import (
        faces_to_grid, subset_model_source, uxgrid_from_corner_bounds)

    lon_b, lat_b = _swath_corner_bounds(obsobj)
    if lon_b is None:
        raise ValueError(
            "_conservative_mod2swath: swath has no corner bounds "
            "(latitude_bounds/longitude_bounds). Conservative regridding "
            f"needs them. Available: {list(obsobj.variables)}."
        )

    olon = np.asarray(obsobj["lon"].values)
    olat = np.asarray(obsobj["lat"].values)
    nx, ny = olon.shape

    # Flatten (x, y, corner) -> (x*y, corner) in C order; the n_face order of
    # the resulting mesh is row-major over (x, y), so reshape back the same way.
    clon = np.asarray(lon_b.values).reshape(nx * ny, -1)
    clat = np.asarray(lat_b.values).reshape(nx * ny, -1)
    swath_grid = uxgrid_from_corner_bounds(clon, clat)

     # model grid file (e.g. SCRIP) or the native MPAS mesh
    grid_file = (modobj.attrs.get("mio_scrip_file")
                 or modobj.attrs.get("mio_mesh_file"))

    src, src_grid = subset_model_source(
        modobj, grid_file,
        float(np.nanmin(clon)), float(np.nanmax(clon)),
        float(np.nanmin(clat)), float(np.nanmax(clat)),
        label="_conservative_mod2swath",
    )
    out = regrid(src, method=method, src_grid=src_grid, target_grid=swath_grid)
    
    result = faces_to_grid(out, (nx, ny), ("x", "y"))
    result.attrs = dict(modobj.attrs)
    return result.assign_coords(
        lon=(("x", "y"), olon),
        lat=(("x", "y"), olat),
    )

def _unstructured_back_to_modgrid(
    concatenated, modobj, method="radius_mean", radius_deg=0.1,):
    
    """Project swath-paired data onto an unstructured model's columns.
    Two modes:

    - ``method`` in (conservative, bilinear, ...) and swath corner bounds
      are present -> true area-weighted mesh-to-mesh regrid via xregrid
      (:func:`_conservative_swath2mod`). **Mass-conserving.**
    - otherwise -> cKDTree within-radius mean: each model column gets the
      mean of all swath pixels within ``radius_deg`` (poor-man's area
      weighting; columns with no pixel in range -> NaN). 0.1 deg ~= ~10 km.
    """
    from melodies_monet.util.regrid_util import regrid

    if method in _UNSTRUCT_CONSERVATIVE:
        lon_b, _ = _swath_corner_bounds(concatenated)
        if lon_b is not None:
            return _conservative_swath2mod_scan([concatenated], modobj, method=method)
        warnings.warn(
            "_unstructured_back_to_modgrid: conservative requested but the "
            "paired swath has no corner bounds; falling back to radius_mean. "
            "(Bounds are carried by _carry_swath_bounds in "
            "regrid_and_apply_weights -- check they survived.)"
        )
    elif method not in _UNSTRUCT_NEAREST:
        # bilinear / patch / anything else hard error with guidance.
        raise _unsupported_method_error(method, "_unstructured_back_to_modgrid")

    # Flatten swath (x, y) -> 1-D "pixel" with longitude/latitude as coords
    # so regrid() sees it as an unstructured source.
    flat = (
        concatenated
        .stack(pixel=("x", "y"))
        .reset_index("pixel", drop=True)
    )
    if "lon" in flat.variables and "longitude" not in flat.variables:
        flat = flat.rename({"lon": "longitude", "lat": "latitude"})
    flat = flat.set_coords(["longitude", "latitude"])

    mlon = np.asarray(modobj["longitude"].values).ravel()
    mlat = np.asarray(modobj["latitude"].values).ravel()
    col_dim = modobj["longitude"].dims[0]

    return regrid(
        flat,
        target={"lon": mlon, "lat": mlat},
        method="radius_mean",
        radius=radius_deg,
        target_dims=(col_dim,),
    )

def tempo_interp_mod2swath(obsobj, modobj, method="conservative", weights=None):
    """Interpolate model to satellite swath/swaths

    Parameters
    ----------
    obsobj : xr.Dataset
        satellite with swath data.
    modobj : xr.Dataset
        model data (with no2 col calculated)
    method : str
        Choose regridding method. Can be "conservative", "conservative_normed",
        "bilinear" or "patch". Check xesmf documentation for details.
    weights : str
        Path to the weightfile. If present, the weights won't be calculated again.

    Returns
    -------
    xr.Dataset
        Regridded model data at swath or swaths. If type is xr.Dataset, a single
        swath is returned. If type is collections.OrderedDict, it returns an
        OrderedDict in which each time represents the reference time of the swath.
    """

    mod_at_swathtime = modobj.interp(time=obsobj.time.mean())

    if mod_at_swathtime.attrs.get("mio_has_unstructured_grid", False):
        # Unstructured model. conservative and true area-weighted mesh-to-mesh
        # via xregrid (mass-conserving). nearest goes via fast cKDTree. bilinear/
        # patch are NOT supported (ESMF needs node-located data).
        if method in _UNSTRUCT_CONSERVATIVE:
            return _conservative_mod2swath(mod_at_swathtime, obsobj, method=method)
        
        if method in _UNSTRUCT_NEAREST:
            return _nearest_mod2swath(mod_at_swathtime, obsobj)
        raise _unsupported_method_error(method, "tempo_interp_mod2swath")
                    
    if weights is None:
        regridder = xe.Regridder(
            mod_at_swathtime,
            obsobj,
            method,
            ignore_degenerate=True,
            unmapped_to_nan=True,
        )
        modswath = regridder(mod_at_swathtime)
    else:
        regridder = xe.Regridder(
            mod_at_swathtime,
            obsobj,
            method,
            ignore_degenerate=True,
            unmapped_to_nan=True,
            filename=weights,
            reuse_weights=True,
        )
        modswath = regridder(mod_at_swathtime)

    return modswath


def _calc_dp(obsobj):
    """Calculate delta pressure in satellite layers

    Parameters
    ----------
    obsobj : xr.Dataset
        satellite observations containing pressure (in Pa)

    Returns
    -------
    xr.DataArray
        Pressure difference in layer
    """

    # REMINDER: pressure is higher at lower vertical levels: dp is positive
    # only if defined as lower - higher.
    dp_vals = (
        obsobj["pressure"].isel(swt_level_stagg=slice(None, -1)).values
        - obsobj["pressure"].isel(swt_level_stagg=slice(1, None)).values
    )
    dp = xr.DataArray(
        data=dp_vals,
        dims=("swt_level", "x", "y"),
        coords={
            "lon": (("x", "y"), obsobj["lon"].values),
            "lat": (("x", "y"), obsobj["lat"].values),
        },
        attrs={
            "units": "Pa",
            "description": "Delta pressure in layer",
            "long_name": "delta_p",
        },
    )
    return dp


@numba.jit(nopython=True)
def _interp_vert(orig, target, data):
    """Performs the numpy interpolation. It is separated from other functions
    for the sake of using the numba jit.

    Parameters:
    -----------
    orig : np.ndarray
        Original grid from which to interpolate. The expected dimensions are (z, x, y),
        in that order. The horizontal and time dimensions are expected to be previously
        interpolated. The original pressure levels should be in decreasing order.
    target : np.ndarray
        Target data with vertical grid information. The expected dimensions are (z, x, y),
        in that order. The target pressure layers should be in decreasing order.
    data : np.ndarray
        Data to be interpolated. It should have the same grid (including vertical) and dimensions
        as orig.

    Returns
    -------
    np.ndarray
        Interpolated data
    """
    assert orig.shape == data.shape, "Grid shape does not match data"
    nz, nx, ny = target.shape
    interp = np.zeros((nz, nx, ny))
    for x in range(nx):
        for y in range(ny):
            interp[:, x, y] = np.flip(
                np.interp(
                    np.flip(target[:, x, y]),
                    np.flip(orig[:, x, y]),
                    np.flip(data[:, x, y]),
                )
            )
    return interp


def calc_altitude_from_thickness(dz_m):
    """Calculate layer altitude above ground

    Parameters
    ----------
    data : xr.DataArray
        DataArray containing dz_m

    Returns
    -------
    Model altitude in satellite space
    """
    altitude_interface = xr.zeros_like(dz_m)
    altitude_interface[{"z": 0}] = dz_m[{"z": 0}]
    for lev in altitude_interface["z"][1:]:
        altitude_interface[{"z": lev}] = altitude_interface[{"z": lev - 1}] + dz_m[{"z": lev}]
    altitude_interface.attrs = {
        "description": "Altitude AGL in m at layer interface",
        "units": "m",
        "long_name": "altitude_agl",
    }
    return altitude_interface


def calc_dz_m_from_altitude(altitude):
    """Calculates dz_m from altitude AGL (in m).

    Parameters
    ----------
    altitude : xr.DataArray
        DataArray containing the layer interface altitude AGL at the
        interface.

    Returns
    -------
    xr.DataArray
        DataArray containing the layer thickness (dz_m) in m.
    """
    dz_m = xr.zeros_like(altitude)
    dz_m[{"z": 0}] = altitude[{"z": 0}]
    dz_m[{"z": slice(1, None)}] = (
        altitude[{"z": slice(1, None)}].values - altitude[{"z": slice(0, -1)}].values
    )
    dz_m.attrs = {
        "description": "Layer thickness in m",
        "units": "m",
        "long_name": "layer_thickness",
    }
    return dz_m


def interp_vertical_mod2swath(obsobj, modobj, variables="NO2_col"):
    """Interpolates model vertical layers to TEMPO vertical layers

    Parameters
    ----------
    modobj : xr.Dataset
        Model data (as provided by MONETIO)
    obsobj : xr.Dataset
        TEMPO data (as provided by MONETIO). Must include pressure.
    variables : str | list[str]
        Variables to interpolate.

    Returns
    -------
    xr.Dataset
        Model data (interpolated to TEMPO vertical layers
    """
    assert np.all(modobj["lon"].fillna(0).values == obsobj["lon"].fillna(0).values)
    assert np.all(modobj["lat"].fillna(0).values == obsobj["lat"].fillna(0).values)

    modsatlayers = xr.Dataset()
    p_mid_tempo = (
        obsobj["pressure"].isel(swt_level_stagg=slice(None, -1)).values
        + obsobj["pressure"].isel(swt_level_stagg=slice(1, None)).values
    ) / 2
    p_orig = modobj["pres_pa_mid"].values
    dimensions = ("z", "x", "y")
    coords = {
        "lon": (("x", "y"), modobj["lon"].values),
        "lat": (("x", "y"), modobj["lat"].values),
    }
    for var in list(variables):
        interpolated = _interp_vert(p_orig, p_mid_tempo, modobj[var].values)
        modsatlayers[var] = xr.DataArray(
            data=interpolated, dims=dimensions, coords=coords, attrs=modobj[var].attrs
        )
    modsatlayers["pres_pa_mid"] = xr.DataArray(
        data=p_mid_tempo,
        dims=dimensions,
        coords=coords,
        attrs=modobj["pres_pa_mid"].attrs,
    )
    _interp_description = "Mid layer pressure interpolated to TEMPO mid swt_layer pressures"
    modsatlayers["pres_pa_mid"].attrs["description"] = _interp_description
    return modsatlayers


def calc_partialcolumn(modobj, var="NO2"):
    """Calculates the partial column of a species from its concentration.

    Parameters
    ----------
    modobj : xr.Dataset
        Model data
    var : str
        variable to calculate the partial column from

    Returns
    -------
    xr.DataArray
        DataArray containing the partial column of the species.
    """
    R = 8.314  # m3 * Pa / K / mol
    NA = 6.022e23
    ppbv2molmol = 1e-9
    m2_to_cm2 = 1e4
    fac_units = ppbv2molmol * NA / m2_to_cm2
    partial_col = (
        modobj[var]
        * modobj["pres_pa_mid"]
        * modobj["dz_m"]
        * fac_units
        / (R * modobj["temperature_k"])
    )
    return partial_col


def apply_weights_mod2tempo_no2_hydrostatic(obsobj, modobj, species="NO2"):
    """Apply the scattering weights and air mass factors according to
    Cooper et. al, 2020, doi: https://doi.org/10.5194/acp-20-7231-2020,
    assuming the hydrostatic equation. It does not require temperature
    nor geometric layer thickness.

    Parameters
    ----------
    obsobj : xr.Dataset
        TEMPO data, including pressure and scattering weights
    modobj : xr.Dataset
        Model data, already interpolated to TEMPO grid

    Returns
    -------
    xr.DataArray
        A xr.DataArray containing the NO2 model data after applying
        the air mass factors and scattering weights
    """
    unit_c = 6.022e23 * 9.8 / 1e4  # NA * g / m2_to_cm2
    dp = _calc_dp(obsobj).rename({"swt_level": "z"})
    ppbv2molmol = 1e-9
    tropopause_pressure = obsobj["tropopause_pressure"]
    scattering_weights = obsobj["scattering_weights"].transpose("swt_level", "x", "y")
    scattering_weights = scattering_weights.rename({"swt_level": "z"})
    scattering_weights = scattering_weights.where(modobj["pres_pa_mid"] >= tropopause_pressure)
    modno2 = modobj[species].where(modobj["pres_pa_mid"] >= tropopause_pressure)
    amf_troposphere = obsobj["amf_troposphere"]
    modno2col_trfmd = (dp * scattering_weights * modno2).sum(dim="z") * unit_c * ppbv2molmol
    modno2col_trfmd = modno2col_trfmd.where(modno2.isel(z=0).notnull())
    modno2col_trfmd = modno2col_trfmd / amf_troposphere
    modno2col_trfmd.attrs = {
        "units": "molecules/cm2",
        "description": "model NO2 tropospheric column after applying TEMPO scattering weights and AMF",
        "history": "Created by MELODIES-MONET, apply_weights_mod2tempo_no2_hydrostatic, TEMPO util",
    }
    return modno2col_trfmd.where(np.isfinite(modno2col_trfmd))


def apply_weights_mod2tempo_no2(obsobj, modobj, species="NO2", column_type="tropospheric"):
    """Apply the scattering weights and air mass factors according to
    Cooper et. al, 2020, doi: https://doi.org/10.5194/acp-20-7231-2020

    Parameters
    ----------
    obsobj : xr.Dataset
        TEMPO data, including pressure and scattering weights
    modobj : xr.Dataset
        Model data, already interpolated to TEMPO grid
    column_type : str
        Whether the "tropospheric" or "total" column should be used for the calculation

    Returns
    -------
    xr.DataArray
        A xr.DataArray containing the NO2 model data after applying
        the air mass factors and scattering weights
    """
    partial_col = modobj[f"{species}_col"]

    scattering_weights = obsobj["scattering_weights"].transpose("swt_level", "x", "y")
    scattering_weights = scattering_weights.rename({"swt_level": "z"})
    if column_type == "tropospheric":
        tropopause_pressure = obsobj["tropopause_pressure"] * 100
        scattering_weights = scattering_weights.where(modobj["pres_pa_mid"] >= tropopause_pressure)
        amf_troposphere = obsobj["amf_troposphere"]
    modno2col_trfmd = (scattering_weights * partial_col).sum(dim="z") / amf_troposphere
    modno2col_trfmd = modno2col_trfmd.where(modobj[f"{species}_col"].isel(z=0).notnull())

    # diagnostic: AK-applied vs RAW tropospheric column. The AK should suppress
    # (ratio < 1)

    # [AK-prof] TEMPO pres_pa_mid (top..bot): [1.e+05 9.e+04 7.e+04 3.e+04 7.e+03 2.e+03 2.e+02 2.e+01] 
    # AK prof pressure being applied backwards?? 
    try:
        _raw = partial_col.where(scattering_weights.notnull()).sum(dim="z")
        _r = (modno2col_trfmd / _raw).where(_raw != 0)
        _swm = float(np.nanmean(scattering_weights.values))
        _amfm = float(np.nanmean(np.asarray(amf_troposphere.values)))
        print(f"[AK] TEMPO {species}: <SW>={_swm:.3f} <AMF_trop>={_amfm:.4f} "
              f"SW/AMF={_swm / _amfm:.2f} | raw_col={float(np.nanmean(_raw.values)):.2e} "
              f"AK-applied={float(np.nanmean(modno2col_trfmd.values)):.2e} "
              f"AK/raw={float(np.nanmean(_r.values)):.2f} (=AMF_mod/AMF_ret)", flush=True)
    except Exception as _e:  # noqa: BLE001
        print(f"[AK] TEMPO {species} ratio diag skipped: {_e!r}", flush=True)
        
    modno2col_trfmd.attrs = {
        "units": "molecules/cm2",
        "description": "model NO2 tropospheric column after applying TEMPO scattering weights and AMF",
        "history": "Created by MELODIES-MONET, apply_weights_mod2tempo_no2, TEMPO util",
    }
    return modno2col_trfmd.where(np.isfinite(modno2col_trfmd))


def apply_weights_mod2tempo_hcho_hydrostatic(obsobj, modobj, species="HCHO"):
    """Apply the scattering weights and air mass factors according to
    Cooper et. al, 2020, doi: https://doi.org/10.5194/acp-20-7231-2020,
    assuming the hydrostatic equation. It does not require temperature
    nor geometric layer thickness.

    Parameters
    ----------
    obsobj : xr.Dataset
        TEMPO data, including pressure and scattering weights
    modobj : xr.Dataset
        Model data, already interpolated to TEMPO grid

    Returns
    -------
    xr.DataArray
        A xr.DataArray containing the NO2 model data after applying
        the air mass factors and scattering weights
    """
    unit_c = 6.022e23 * 9.8 / 1e4  # NA * g / m2_to_cm2
    dp = _calc_dp(obsobj).rename({"swt_level": "z"})
    ppbv2molmol = 1e-9
    scattering_weights = obsobj["scattering_weights"].transpose("swt_level", "x", "y")
    scattering_weights = scattering_weights.rename({"swt_level": "z"})
    modhcho = modobj[species]
    amf = obsobj["amf"]
    modhcho_col = (dp * scattering_weights * modhcho).sum(dim="z") * unit_c * ppbv2molmol
    modhcho_col = modhcho_col / amf
    modhcho_col.attrs = {
        "units": "molecules/cm2",
        "description": "model HCHO column after applying TEMPO scattering weights and AMF",
        "history": "Created by MELODIES-MONET, apply_weights_mod2tempo_hcho_hydrostatic, TEMPO util",
    }
    return modhcho_col.where(np.isfinite(modhcho_col))


def apply_weights_mod2tempo_hcho(obsobj, modobj, species="HCHO"):
    """Apply the scattering weights and air mass factors according to
    Cooper et. al, 2020, doi: https://doi.org/10.5194/acp-20-7231-2020

    Parameters
    ----------
    obsobj : xr.Dataset
        TEMPO data, including pressure and scattering weights
    modobj : xr.Dataset
        Model data, already interpolated to TEMPO grid

    Returns
    -------
    xr.DataArray
        A xr.DataArray containing the NO2 model data after applying
        the air mass factors and scattering weights
    """
    partial_col = modobj[f"{species}_col"]

    scattering_weights = obsobj["scattering_weights"].transpose("swt_level", "x", "y")
    scattering_weights = scattering_weights.rename({"swt_level": "z"})
    amf = obsobj["amf"]
    modhcho_col = (scattering_weights * partial_col).sum(dim="z") / amf
    modhcho_col = modhcho_col.where(modobj[f"{species}_col"].isel(z=0).notnull())
    modhcho_col.attrs = {
        "units": "molecules/cm2",
        "description": "model HCHO column after applying TEMPO scattering weights and AMF",
        "history": "Created by MELODIES-MONET, apply_weights_mod2tempo_hcho, TEMPO util",
    }
    return modhcho_col.where(np.isfinite(modhcho_col))


def is_nonpairable(obsobj, k, modobj):
    """Discards inplace granules from obsobj that do not match modobj's
    domain, or granules that are all NaN. If the domain is small,
    it can considerably speed up the regridding process.

    Parameters
    ----------
    obsobj : dict[str, xr.Dataset]
        tempo data
    modobj : xr.Dataset
        model data

    Return
    ------
    collections.OrderedDict[str, xr.Dataset]
    """
    
    # if obsobj[k]["lon"].max() < modobj["longitude"].min():
    
    # Normalize longitudes to the same [-180, 180] convention before
    # comparing. CESM-SE (and several other models) store lon in 0..360, while
    # TEMPO L2 reports in -180..180; a CONUS model and CONUS
    # swath look non-overlapping (e.g. model min=200 vs obs max=-60) and every
    # granule is discarded.
    
    def _to_neg180(a):
        import numpy as np
        return ((np.asarray(a) + 180.0) % 360.0) - 180.0

    mlon = _to_neg180(modobj["longitude"].values)
    olon = _to_neg180(obsobj[k]["lon"].values)
            
    if olon.max() < mlon.min():
        return True
    elif olon.min() > mlon.max():
        return True

    #elif obsobj[k]["lon"].min() > modobj["longitude"].max():
    elif olon.min() > mlon.max():
        return True
    elif obsobj[k]["lat"].max() < modobj["latitude"].min():
        return True
    elif obsobj[k]["lat"].min() > modobj["latitude"].max():
        return True
    return False


def _regrid_and_apply_weights(
    obsobj, modobj, method="conservative", weights=None, species=["NO2"], tempo_sp="NO2",
    crop_extent=None,
):
    """Does the complete process of regridding and
    applying scattering weights. Assumes that obsobj is a Dataset

    Parameters
    ----------
    obsobj : xr.Dataset
        TEMPO observations
    modobj : xr.Dataset
        Model data
    method : str
        Choose regridding method. Can be "conservative", "conservative_normed",
        "bilinear" or "patch". Check xesmf documentation for details.
    weights : str
        Path to the weightfile. If present, the weights won't be calculated again.
    tempo_sp: str
        NO2 or HCHO, to apply the correct Air Mass Factors and scattering weights.

    Returns
    -------
    xr.DataArray
        Model data regridded to the TEMPO grid,
        with the averaging kernel.
    """
    if tempo_sp == "NO2":
        apply_weights = apply_weights_mod2tempo_no2
        apply_weights_hydrostatic = apply_weights_mod2tempo_no2_hydrostatic
    else:
        assert tempo_sp == "HCHO", "TEMPO species must be HCHO or NO2."
        apply_weights = apply_weights_mod2tempo_hcho
        apply_weights_hydrostatic = apply_weights_mod2tempo_hcho_hydrostatic
    if method == "conservative" and not modobj.attrs.get(
        "mio_has_unstructured_grid", False
    ):
        if "lat_b" not in modobj:
            calc_grid_corners(modobj)
        if "lat_b" not in obsobj:
            calc_grid_corners(obsobj, lat="lat", lon="lon")
    modobj_hs = tempo_interp_mod2swath(obsobj, modobj, method=method, weights=weights)
    if "dz_m" in modobj.keys():
        modobj_hs["altitude"] = calc_altitude_from_thickness(modobj_hs["dz_m"])
        modobj_swath = interp_vertical_mod2swath(
            obsobj, modobj_hs, [f"{species[0]}", "altitude", "temperature_k"]
        )
        modobj_swath["dz_m"] = calc_dz_m_from_altitude(modobj_swath["altitude"])
        modobj_swath[f"{species[0]}_col"] = calc_partialcolumn(modobj_swath, var=species[0])

        # more debugging

        # [AK-cons] column conservation across the vertical interp: the
        # tropospheric column on NATIVE model levels vs after interp to TEMPO
        # levels must match (~1); a big departure = the interp misplaces mass,
        # while ~1 means any AK/raw > 1 is profile-shape physics.
        
        try:
            _sp = species[0]
            # bug in xeregrid? 
            # # pre-regrid (straight from the reader) vs post-regrid surface values.
            # # If post = 2x pre, the horizontal regrid double-counts; if T is also
            # # ~2x (~580 K), EVERY field is doubled (unnormalized weights).
            # _p0v = float(np.nanmax(np.asarray(modobj["pres_pa_mid"].values)))
            # _p1v = float(np.nanmax(np.asarray(modobj_hs["pres_pa_mid"].values)))
            # _t1v = float(np.nanmean(np.asarray(modobj_hs["temperature_k"].isel(z=0).values)))
            # print(f"[AK-rgd] pres max pre-regrid={_p0v:.0f} Pa | post-regrid={_p1v:.0f} Pa "
            #       f"(ratio {_p1v / _p0v:.2f}) | post-regrid sfc T={_t1v:.0f} K (expect ~290)",
            #       flush=True)
            # _pn = np.nanmean(np.asarray(modobj_hs["pres_pa_mid"].values), axis=(1, 2))
            # _xn = np.nanmean(np.asarray(modobj_hs[_sp].values), axis=(1, 2))
            # _pi = np.nanmean(np.asarray(modobj_swath["pres_pa_mid"].values), axis=(1, 2))
            # _xi = np.nanmean(np.asarray(modobj_swath[_sp].values), axis=(1, 2))
            # _k = max(len(_pn) // 8, 1)
            # _ki = max(len(_pi) // 8, 1)
            # print(f"[AK-nat] native pres : {np.array2string(_pn[::_k], precision=0)}", flush=True)
            # print(f"[AK-nat] native {_sp} : {np.array2string(_xn[::_k], precision=3)}", flush=True)
            # print(f"[AK-int] interp pres : {np.array2string(_pi[::_ki], precision=0)}", flush=True)
            # print(f"[AK-int] interp {_sp} : {np.array2string(_xi[::_ki], precision=3)}", flush=True)
            
            # # [AK-cons] column conservation across the vertical interp: the
            # # tropospheric column on NATIVE model levels vs after interp to
            # # TEMPO levels must match (~1). A big departure = the interp is
            # # misplacing mass; ~1 = any AK/raw > 1 is profile-shape physics.
            _pc_nat = calc_partialcolumn(modobj_hs, var=_sp)
            if "tropopause_pressure" in obsobj:
                _tp = obsobj["tropopause_pressure"] * 100
            else:
                _tp = 0.0 * modobj_hs["pres_pa_mid"].isel(z=0)  # full column
            _cn = _pc_nat.where(modobj_hs["pres_pa_mid"] >= _tp).sum("z")
            _ci = (modobj_swath[f"{_sp}_col"]
                   .where(modobj_swath["pres_pa_mid"] >= _tp).sum("z"))
            _cr = (_ci / _cn).where(_cn != 0)
            print(f"[AK-cons] {_sp} trop col native={float(np.nanmean(_cn.values)):.2e} "
                  f"interp={float(np.nanmean(_ci.values)):.2e} "
                  f"interp/native mean={float(np.nanmean(_cr.values)):.2f} "
                  f"median={float(np.nanmedian(_cr.values)):.2f} (expect ~1)", flush=True)
            
        except Exception as _e:  # noqa: BLE001
            print(f"[AK-cons] skipped: {_e!r}", flush=True)
        da_out = apply_weights(obsobj, modobj_swath, species=f"{species[0]}")
    else:
        warnings.warn(
            "There is no dz_m variable, and the partial column"
            + "cannot be directly calculated. Assuming hydrostatic equation."
        )
        modobj_swath = interp_vertical_mod2swath(obsobj, modobj_hs, species)
        da_out = apply_weights_hydrostatic(obsobj, modobj_swath, species=species[0])
    return da_out.where(np.isfinite(da_out))


def regrid_and_apply_weights(
    obsobj,
    modobj,
    pair=True,
    verbose=True,
    method="conservative",
    weights=None,
    species=["NO2"],
    tempo_sp="NO2",
    crop_extent=None,
):
    """Does the complete process of regridding
    and applying scattering weights.

    Parameters
    ----------
    obsobj : xr.Dataset | collections.OrderedDict
        TEMPO observations
    modobj : xr.Dataset
        Model output
    pair : boolean
        If True, returns paired data.
    verbose : boolean
        If True, let's the user know when each timestamp is being regridded.
        Only has an effect if the input is an OrderedDict
    method : str
        Choose regridding method. Can be "conservative", "conservative_normed",
        "bilinear" or "patch". Check xesmf documentation for details.
    weights : None | str
        If present, a weightfile (as in "weights") is applied
    discard_useless: boolean
        If True, satellite granules that don't match the model domain are not used.
    tempo_sp: str
        NO2 for the NO2 product, HCHO for the HCHO product

    Returns
    -------
    xr.Dataset | collections.OrderedDict
        Model with regridded data. If obsobj is of type collections.OrderedDict,
        an OrderedDict is returned.
    """

    if tempo_sp == "NO2":
        sat_species_name = "vertical_column_troposphere"
    else:
        assert tempo_sp == "HCHO", "TEMPO species must be HCHO or NO2."
        sat_species_name = "vertical_column"

    def _obs_pair_vars(src):
        """Obs column + any per-pixel retrieval uncertainty/precision that goes
        with it. Carrying the error forward lets 'series' inverse-variance
        weight (instead of silently falling back to area) and lets 'swath'
        preserve per-pixel uncertainty. 
        """
        names = [sat_species_name]
        for v in src.variables:
            if v == sat_species_name:
                continue
            if v.startswith(sat_species_name) and v.endswith(("_uncertainty", "_precision")):
                names.append(v)
        return src[names]
        
    # crop granule to region of interest
    from melodies_monet.util.sat_l2_swath_utility import _crop_swath_to_extent
    
    if isinstance(obsobj, xr.Dataset):
        if crop_extent is not None:
            obsobj = _crop_swath_to_extent(obsobj, crop_extent)
            if obsobj is None:
                return None
                
        regridded = _regrid_and_apply_weights(
            obsobj, modobj, method=method, weights=weights, species=species, tempo_sp=tempo_sp
        )
        output = regridded.to_dataset(name=species[0])
        output.attrs["reference_time_string"] = obsobj.attrs["reference_time_string"]
        output.attrs["final_time_string"] = obsobj["time"][-1].values.astype(str)
        output.attrs["scan_num"] = obsobj.attrs["scan_num"]
        output.attrs["granule_number"] = obsobj.attrs["granule_number"]
        if pair:
            output = xr.merge([output, _obs_pair_vars(obsobj)])
        if "lat" in output.variables:
            output = output.rename({"lat": "latitude", "lon": "longitude"})
        output = _carry_swath_bounds(output, obsobj)
        return output
    if isinstance(obsobj, dict):
        output_multiple = {}
        for ref_time in obsobj.keys():
            if is_nonpairable(obsobj, ref_time, modobj):
                warnings.warn(f"{ref_time} granule domain has no overlap with model. Discarding.")
                continue
                
            granule = obsobj[ref_time]
            if crop_extent is not None:
                _pre = dict(granule.sizes)
                granule = _crop_swath_to_extent(granule, crop_extent)
                if granule is None:
                    if verbose:
                        print(f"  {ref_time}: no overlap with crop_extent="
                              f"{crop_extent}; granule discarded", flush=True)
                    continue
                if verbose:
                    print(f"  {ref_time}: cropped {_pre} -> "
                          f"{dict(granule.sizes)}", flush=True)
            elif verbose:
                print(f"  {ref_time}: crop_extent is None; granule NOT cropped",
                      flush=True)
            if verbose:
                print(f"Regridding {ref_time} and applying AMF and weights")
            output_multiple[ref_time] = _regrid_and_apply_weights(
                granule,
                modobj,
                method=method,
                weights=weights,
                species=species,
                tempo_sp=tempo_sp,
            ).to_dataset(name=species[0])
            output_multiple[ref_time].attrs["reference_time_string"] = ref_time
            output_multiple[ref_time].attrs["scan_num"] = obsobj[ref_time].attrs["scan_num"]
            output_multiple[ref_time].attrs["granule_number"] = obsobj[ref_time].attrs[
                "granule_number"
            ]
            output_multiple[ref_time].attrs["final_time_string"] = obsobj[ref_time]["time"][
                -1
            ].values.astype(str)
            if pair:
                output_multiple[ref_time] = xr.merge(
                    [
                        output_multiple[ref_time],
                        _obs_pair_vars(granule),  # obs col + uncertainty, cropped 
                    ]
                )
            if "lat" in output_multiple[ref_time].variables:
                output_multiple[ref_time] = output_multiple[ref_time].rename(
                    {"lat": "latitude", "lon": "longitude"}
                )

            # better memory handling 
            output_multiple[ref_time] = _carry_swath_bounds(
                output_multiple[ref_time], granule
            )
            
            # better memory handling 
            output_multiple[ref_time] = output_multiple[ref_time].load()
            gc.collect()

            import psutil, os
            proc = psutil.Process(os.getpid())
            print(f"  rss={proc.memory_info().rss / 1e9:.2f} GB")

        return output_multiple
    raise TypeError("Obsobj must be xr.Dataset or dict")

def _conservative_swath2mod_scan(granules, modobj, method="conservative"):
    """Conservative swath to unstructured-mesh regrid, accumulated PER GRANULE.

    Regridding a whole scan at once builds a (huge swath mesh x full model
    mesh) weight matrix that OOMs
    
    Regrid each granule separately against a model-mesh SUBSET (the granule's bbox) and coverage-accumulate
    the NON-normalized conservative outputs onto the full mesh::

        combined_mean = sum_g(raw_g) / sum_g(cov_g)

    where ``raw_g`` is the (destarea) regrid of the field for granule g and
    ``cov_g`` is the regrid of ones (fraction of each model cell covered by
    that granule). Both source (one granule) and target (a strip) stay small,
    so it scales to the full ne0CONUS mesh without OOM. Uncovered cells
    (sum cov ~ 0) become NaN.
    
    """
    import uxarray as ux
    from melodies_monet.util.regrid_util import regrid
    from melodies_monet.util.uxarray_util import (
        flatten_to_faces, open_uxgrid, subset_mesh_to_bbox,
        uxgrid_from_corner_bounds)

    _grid_file = (modobj.attrs.get("mio_scrip_file")
                  or modobj.attrs.get("mio_mesh_file"))
    model_grid = open_uxgrid(_grid_file)
    n_col_full = int(model_grid.n_face)

    num = {}                                       # var full-mesh numerator
    attrs_by_var = {}
    den = np.zeros(n_col_full, dtype="float64")    # coverage accumulator

    for g in granules:
        lon_b, lat_b = _swath_corner_bounds(g)
        if lon_b is None:
            continue
        olon = np.asarray(g["longitude"].values if "longitude" in g.variables
                          else g["lon"].values)
        nx, ny = olon.shape
        clon = np.asarray(lon_b.values).reshape(nx * ny, -1)
        clat = np.asarray(lat_b.values).reshape(nx * ny, -1)
        swath_grid = uxgrid_from_corner_bounds(clon, clat)

        flat = flatten_to_faces(
            g, ("x", "y"),
            drop=(lon_b.name, lat_b.name, "longitude", "latitude", "lon", "lat"))
        
        flat = flat.assign(
            _mm_cov=(("n_face",), np.ones(flat.sizes["n_face"], dtype="float64")))
        src_uxds = ux.UxDataset(flat, uxgrid=swath_grid)

        # target subset = model faces in this granule's bbox
        _keep, tgt = subset_mesh_to_bbox(
            _grid_file,
            float(np.nanmin(clon)), float(np.nanmax(clon)),
            float(np.nanmin(clat)), float(np.nanmax(clat)))
        if _keep.size == 0:
            continue
        if tgt is None:
            _keep, tgt = np.arange(n_col_full), model_grid
        out_sub = regrid(src_uxds, method=method, target_grid=tgt)
        _sfd = next((d for d in out_sub.dims if out_sub.sizes[d] == _keep.size), None)
        if _sfd is None:
            continue
        den[_keep] += np.nan_to_num(
            np.asarray(out_sub["_mm_cov"].transpose(_sfd).values))
        for v in out_sub.data_vars:
            if v == "_mm_cov":
                continue
            da = out_sub[v].transpose(..., _sfd)
            arr = np.nan_to_num(np.asarray(da.values))
            if v not in num:
                num[v] = np.zeros(arr.shape[:-1] + (n_col_full,), dtype="float64")
                attrs_by_var[v] = dict(out_sub[v].attrs)
            num[v][..., _keep] += arr
        print(f"_conservative_swath2mod_scan: granule nx*ny={nx * ny} | "
              f"target {n_col_full} -> {_keep.size} faces", flush=True)
        del out_sub, src_uxds
        gc.collect()

    col_dim = modobj["longitude"].dims[0]
    den_ok = den > 1e-6
    out = xr.Dataset()
    for v, arr in num.items():
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = arr / den
        mean = np.where(den_ok, mean, np.nan)
        extra = tuple(f"_stack{i}" for i in range(arr.ndim - 1))
        out[v] = xr.DataArray(mean, dims=extra + (col_dim,), attrs=attrs_by_var[v])
    out = out.assign_coords({
        "longitude": (col_dim, np.asarray(modobj["longitude"].values).ravel()),
        "latitude": (col_dim, np.asarray(modobj["latitude"].values).ravel()),
    })
    return out

def back_to_modgrid(
    paireddict,
    modobj,
    keys_to_merge="all",
    add_time=True,
    to_netcdf=False,
    path="Regridded_object_XYZ.nc",
    method="bilinear",
    grid_path=None,
    regrid_target="model",
    obs_grid_res=0.1,
    obs_grid_units="deg",
    obs_grid_extent=None,
):
    """Grids object in sat-space to modgrid. Designed to grid back to modgrid after applying
    the scattering weights and air mass factors. It is designed for a single scan.

    Parameters
    ----------
    paireddict : collections.OrderedDict[str, xr.Dataset]
        An OrderedDict with time_reference strings as keys.
    modobj : xr.Dataset
        A modobj including the modgrid.
    keys_to_merge : str | list[str]
        If 'all', all keys are assumed to be part of the same scan and merged.
        Else, only the keys provided are merged.
    add_time : bool
        If True, add reference time as a coordinate for the scan.
        Can be useful to concatenate later if multiple scans are required.
    to_netcdf : bool
        If True, save a netcdf with the paired data
    path : str
        The base name to save the files if to_netcdf is True. XX will be replaced
        with the scan number and the reference time. If to_netcdf is False, this will
        be ignored.
    method : str
        Method of regridding used by xESMF
    grid_path : str
        If None, defaults to the model grid. Otherwise, the grid in path is used.
        If the method is conservative, lat_b and lon_b are required.

    Returns
    -------
    xr.Dataset
        Dataset with obj2grid regridded to modobj.
    """
    if keys_to_merge == "all":
        ordered_keys = sorted(list(paireddict.keys()))
    else:
        ordered_keys = sorted(list(keys_to_merge))

    def _yidx(ds):
        # Give the across-track dim a positional index so sibling granules of a
        # scan cropped to DIFFERENT across-track windows can be outer-joined at
        # concat instead of raising an unindexed-dim size mismatch
        # Every pixel keeps its own 2-D lon/lat, so per-pixel regrids are unaffected by the
        # positional alignment; padded cells are NaN and ignored downstream.
        if "y" in ds.dims and "y" not in ds.indexes:
            return ds.assign_coords(y=("y", np.arange(ds.sizes["y"])))
        return ds

    concatenated = _yidx(paireddict[ordered_keys[0]])
    
    scan_num = concatenated.attrs["scan_num"]
    granules = [concatenated.attrs["granule_number"]]
    ref_times = [concatenated.attrs["reference_time_string"][:-1]]  # Remove unneeded Z
    if len(paireddict) > 1:
        for k in ordered_keys[1:]:
            ds_to_add = paireddict[k]
            if ds_to_add.attrs["scan_num"] != scan_num:
                raise ValueError(
                    "back_to_modgrid is prepared to work with data of a single scan. "
                    + f"However, {ordered_keys[0]} is from scan {scan_num} and "
                    + f"{k} if from scan {ds_to_add.attrs['scan_num']}."
                )
            concatenated = xr.concat(
                [concatenated, _yidx(paireddict[k])], dim="x", join="outer")
            granules.append(paireddict[k].attrs["granule_number"])
            ref_times.append(paireddict[k].attrs["reference_time_string"][:-1])

    end_time = np.array(
        paireddict[ordered_keys[-1]].attrs["final_time_string"], dtype="datetime64[ns]"
    )

    # enables regridding to obs or model space 
    targets = [regrid_target] if isinstance(regrid_target, str) else list(regrid_target)
    
    from melodies_monet.util.sat_l2_swath_utility import _swath2latlon, _model_lonlat_extent, _swath2series
    _extent = tuple(obs_grid_extent) if obs_grid_extent else _model_lonlat_extent(modobj)
    
    scan_num = concatenated.attrs["scan_num"]

    results = {}
    for _tgt in targets:
        if _tgt == "obs":
            # _dvars = [
            #     v for v in concatenated.data_vars
            #     if v not in ("lon", "lat", "longitude", "latitude")
            #     and "bounds" not in v
            # ]
            # out_regridded = _swath2latlon(
            #     concatenated, _dvars, obs_grid_res, _extent, units=obs_grid_units, method=method,
            # )

            # Regrid EACH granule onto the fixed lat/lon grid and coverage-mean,
            # rather than the outer-join-padded scan concat. Keeps each source
            # small AND avoids NaN-padded corner bounds that would silently knock
            # conservative down to radius_mean inside _swath2latlon.
            _acc_sum, _acc_cnt = None, None
            for _gk in ordered_keys:
                _g = paireddict[_gk]
                _gdv = [
                    v for v in _g.data_vars
                    if v not in ("lon", "lat", "longitude", "latitude")
                    and "bounds" not in v
                ]
                _one = _swath2latlon(
                    _g, _gdv, obs_grid_res, _extent,
                    units=obs_grid_units, method=method,
                )
                _cnt = _one.notnull().astype("float64")
                _one0 = _one.fillna(0.0)
                if _acc_sum is None:
                    _acc_sum, _acc_cnt = _one0, _cnt
                else:
                    _acc_sum = _acc_sum + _one0
                    _acc_cnt = _acc_cnt + _cnt
            out_regridded = (_acc_sum / _acc_cnt).where(_acc_cnt > 0)

        elif _tgt =="series": # time domain vector 
            _errs = [v for v in concatenated.data_vars
                     if v.endswith(("_uncertainty", "_precision"))]
            _dvars = [
                v for v in concatenated.data_vars
                if v not in ("lon", "lat", "longitude", "latitude")
                and "bounds" not in v and v not in _errs
            ]
            _ssw, _obs_var = concatenated, None
            if _errs:
                _ssw = concatenated.assign(_obs_err=concatenated[_errs[0]])
                _base = _errs[0].rsplit("_", 1)[0]
                _obs_var = _base if _base in _dvars else None
            out_regridded = _swath2series(_ssw, _dvars, obs_var=_obs_var)
        elif grid_path is not None:
            grid = xr.open_dataset(grid_path)
            regridder = xe.Regridder(concatenated, grid, method=method, unmapped_to_nan=True)
            out_regridded = regridder(concatenated)
        elif modobj.attrs.get("mio_has_unstructured_grid", False):
            # Unstructured target bypass xESMF (memory issues on 1-D ncol targets)
            # Conservative goes mesh-to-mesh via xregrid; else cKDTree radius-mean.
            # For conservative, regrid PER GRANULE and coverage-accumulate onto
            # the full mesh -- concatenating the whole scan first builds a
            # (huge swath x full mesh) weight matrix that OOMs.
            if (method in _UNSTRUCT_CONSERVATIVE
                    and _swath_corner_bounds(concatenated)[0] is not None):
                out_regridded = _conservative_swath2mod_scan(
                    [paireddict[k] for k in ordered_keys], modobj, method=method
                )
            else:
                out_regridded = _unstructured_back_to_modgrid(
                    concatenated, modobj, method=method
                )
        else:
            regridder = xe.Regridder(concatenated, modobj, method=method, unmapped_to_nan=True)
            out_regridded = regridder(concatenated)

        # shared post-processing 
        for v in out_regridded.variables:
            if v in concatenated.variables:
                out_regridded[v].attrs = concatenated[v].attrs
        out_regridded.attrs["reference_time_string"] = ref_times
        out_regridded.attrs["granules"] = np.array(granules)
        out_regridded.attrs["scan_num"] = scan_num
        out_regridded = out_regridded.where(np.isfinite(out_regridded))
        if add_time:
            time = [np.array(ref_times[0], dtype="datetime64[ns]")]
            da_time = xr.DataArray(
                name="time",
                data=time,
                dims=["time"],
                attrs={"description": "Reference start time of first selected granule in scan."},
                coords={"time": (("time",), time)},
            )
            out_regridded = out_regridded.expand_dims(time=da_time)
            out_regridded["end_time"] = (("time",), [end_time])
            out_regridded["end_time"].attrs = {
                "description": "time at which the last swath of the scan starts"
            }
        if to_netcdf and _tgt == "model":
            if "XYZ" in path:
                out_regridded.to_netcdf(
                    path.replace(
                        "XYZ",
                        f"S{scan_num:03d}_{out_regridded['time'].values.astype(str)[0][0:19]}",
                    )
                )
            else:
                out_regridded.to_netcdf(path)
        results[_tgt] = out_regridded
    return results

def _swath_pixels_from_paireddict(paireddict):
    """Native-swath output: every valid TEMPO pixel across all
    granules as a 1-D ``obs`` vector, carrying ``longitude``/``latitude``/
    ``time`` as per-pixel coords. 
    
    No regridding onto any imposed grid -- the
    most faithful product. All-NaN pixels dropped. 

    Only variables that live purely on the swath spatial dims (x, y) are
    kept, so vertical-level fields and any stray non-spatial dims are skipped.
    """
    pieces = []
    for k in sorted(paireddict):
        g = paireddict[k]
        sp = [d for d in ("x", "y") if d in g.dims]
        if not sp:
            continue
        t = np.array(g.attrs["reference_time_string"][:-1], dtype="datetime64[ns]")
        dvars = [
            v for v in g.data_vars
            if "bounds" not in v
            and v not in ("lon", "lat", "longitude", "latitude")
            and set(g[v].dims) <= set(sp)
        ]
        if not dvars:
            continue
        flat = g[dvars].stack(obs=sp).reset_index("obs", drop=True)
        lon = np.asarray(g["lon"].stack(obs=sp).reset_index("obs", drop=True).values)
        lat = np.asarray(g["lat"].stack(obs=sp).reset_index("obs", drop=True).values)
        n = flat.sizes["obs"]
        flat = flat.assign_coords(
            longitude=("obs", lon), latitude=("obs", lat),
            time=("obs", np.full(n, t)),
        )
        # carry the pixel CORNER bounds (obs, corner) so the footprint
        # polygons can be drawn without resampling
        for _bn in ("longitude_bounds", "latitude_bounds"):
            if _bn in g.variables:
                _bb = g[_bn].stack(obs=sp).reset_index("obs", drop=True)
                _cd = next((d for d in _bb.dims if d != "obs"), None)
                if _cd is not None:
                    flat[_bn] = _bb.transpose("obs", _cd)
                    
        pieces.append(flat)
    if not pieces:
        return xr.Dataset()
    out = xr.concat(pieces, dim="obs")
    # drop pixels that are NaN in every data var (outside cropped granule, masked)
    _sci = [v for v in out.data_vars if "bounds" not in v]
    keep = np.zeros(out.sizes["obs"], dtype=bool)
    for v in _sci:
        keep |= np.isfinite(np.asarray(out[v].values))
    out = out.isel(obs=np.where(keep)[0])
    print(f"_swath_pixels_from_paireddict: {out.sizes.get('obs', 0)} valid pixels "
          f"from {len(paireddict)} granules; vars={list(out.data_vars)}", flush=True)
    return out

def back_to_modgrid_multiscan(
    paireddict,
    modobj,
    to_netcdf=False,
    path="Regridded_object_XYZ.nc",
    method="bilinear",
    grid_path=None,
    regrid_target="model",
    obs_grid_res=0.1,
    obs_grid_units="deg",
    obs_grid_extent=None,
):
    """Grids object in sat-space to modgrid. Designed to grid back to modgrid after applying
    the scattering weights and air mass factors. It is designed for multiple scans, and uses
    back_to_modgrid under the hood. Generally, back_to_modgrid should only be used if
    you can ensure that you are reading only one scan at a time.

    Parameters
    ----------
    paireddict : collections.OrderedDict[str, xr.Dataset]
        An OrderedDict with time_reference strings as keys.
    modobj : xr.Dataset
        A modobj including the modgrid.
    to_netcdf : bool
        If True, save a netcdf with the paired data
    path : str
        The base name to save the files if to_netcdf is True. XX will be replaced
        with the first and list times. If to_netcdf is False, this will
        be ignored.
    method : str
        Method of regridding used by xESMF
    regrid_target: defaults to regridding unstructured grids to model space . If user 
                   specifies [model, obs] in yaml key argument, then regridding of obs will occur to model and obs space 
    
    Returns
    -------
    xr.Dataset
        Dataset with obj2grid regridded to modobj.
    """

    targets = [regrid_target] if isinstance(regrid_target, str) else list(regrid_target)
    out_by = {t: xr.Dataset() for t in targets}

    ordered_keys = sorted(list(paireddict.keys()))
    if not ordered_keys:
        raise ValueError(
            "back_to_modgrid_multiscan received an empty paireddict. "
            "Every granule was discarded as 'no overlap with model' upstream. "
            "Check is_nonpairable() and the model/obs longitude conventions "
            "(model in 0..360 vs obs in -180..180 is a common cause).")

    # "swath" (native pixels) is handled separately -- it is NOT gridded,
    # so it never enters the back_to_modgrid target loop (which would try to
    # build a regridder for it). The remaining targets grid as before.
    grid_targets = [t for t in targets if t != "swath"]

    if grid_targets:
        _kw = dict(
            method=method, grid_path=grid_path, regrid_target=grid_targets,
            obs_grid_res=obs_grid_res, obs_grid_units=obs_grid_units,
            obs_grid_extent=obs_grid_extent,
        )

        scan_num = paireddict[ordered_keys[0]].attrs["scan_num"]
        keys_in_scan = [ordered_keys[0]]
        if len(ordered_keys) > 1:
            for k in ordered_keys[1:]:
                if paireddict[k].attrs["scan_num"] == scan_num:
                    keys_in_scan.append(k)
                else:
                    scan_dict = back_to_modgrid(paireddict, modobj, keys_in_scan, **_kw)
                    for t in grid_targets:
                        out_by[t] = xr.merge([out_by[t], scan_dict[t]])
                    scan_num = paireddict[k].attrs["scan_num"]
                    keys_in_scan = [k]
        scan_dict = back_to_modgrid(paireddict, modobj, keys_in_scan, add_time=True, **_kw)
        for t in grid_targets:
            out_by[t] = xr.merge([out_by[t], scan_dict[t]])

    if "swath" in targets:
        out_by["swath"] = _swath_pixels_from_paireddict(paireddict)

    if to_netcdf:
        for t in targets:
            _p = path if t == "model" else path.replace(".nc", f"_{t}.nc")
            if "XYZ" in _p:
                first_time = out_by[t]["time"][0].values.astype(str)[0:19]
                last_time = out_by[t]["time"][-1].values.astype(str)[0:19]
                out_by[t].to_netcdf(_p.replace("XYZ", f"{first_time}_{last_time}"))
            else:
                out_by[t].to_netcdf(_p)

    return out_by


def read_paired_gridded_tempo_model(path):
    """Reads in paired gridded tempo and model data

    Parameters
    ----------
    path : str or globobject

    Returns
    -------
    xr.Dataset
        combined dataset with paired tempo and gridded data.
    """

    pathlist = sorted(glob.glob(path))
    return xr.open_mfdataset(pathlist)


def read_paired_swath(path):
    """Read in paired swath data

    Parameters
    ----------
    path : str or globobject

    Returns
    -------
    collections.OrderedDict[str, xr.Dataset]
        OrderedDict with datasets containing every swath
    """
    pathlist = sorted(glob.glob(path))
    swaths = collections.OrderedDict()
    for f in pathlist:
        ds = xr.open_dataset(f)
        time = ds.attrs["reference_time_string"]
        swaths[time] = ds
    return swaths


def save_paired_swath(moddict, path="Paired_swath_XYZ.nc"):
    """Saves each swath individually

    Parameters
    ----------
    moddict : collections.OrderedDict[str, xr.Dataset]
        Ordered dict containing all of the paired model and swath.
    path : str
        Path to save the swath. If XYZ is present, it will be replaced
        by the date, number of scan and number of granule.

    Returns
    -------
    None
    """
    for i, k in enumerate(moddict.keys()):
        scan_num = moddict[k].attrs["scan_num"]
        gran_num = moddict[k].attrs["granule_number"]
        if isinstance(path, str):
            if "XYZ" in path:
                pathout = path.replace("XYZ", f"{k[:-1]}_S{scan_num:03d}G{gran_num:03d}")
            else:
                pathout = path.replace(".nc", f"{i}.nc")
        else:
            pathout = f"{i}_paired.nc"
        moddict[k].to_netcdf(pathout)


def select_by_keys(data_names, period="per_scan"):
    """Selects data containing the same scan or day. It does
    so by file name, so it is important that the standard naming
    is not altered.

    Parameters
    ----------
    data_names : list[str]
        list containing the names of the files to select from.
    selection_criteria : str
        str with the selection criteria. If "all", the function does nothing.
        If "per_scan", it will return a list of lists[str], with each inner
        list containing one complete scan.
        If "per_day", it will return a list of lists[str], with each inner
        list containing one complete day.

    Returns
    -------
    list[str] | list[list[str]]
        A sorted list of strings if "all" is provided, a list[list[str]]
        if a different criteria is provided.
    """
    import re

    date_names_sorted = sorted(data_names)

    if period not in ("all", "per_day", "per_scan"):
        warnings.warn(
            "Could not understand looping recommendations for processing files."
            ' Not in "all", "per_day" nor "per_scan". Adopting "all".'
        )
        period = "all"

    if period != "all":
        days = sorted({re.search(r"((\d{8}))T(\d{6})", s).group(1) for s in date_names_sorted})
        subgroups = []
        for day in days:
            subgroups.append([d for d in date_names_sorted if day in d])

        if period == "per_day":
            return subgroups

        if period == "per_scan":
            scans = []
            for sg in subgroups:
                scan = sorted({re.search(r"(S(\d{3}))G(\d{2})", s).group(1) for s in sg})
                for s in scan:
                    scans.append([f for f in sg if s in f])
            return scans

    # This only should happen if period == "all"
    return [date_names_sorted]
