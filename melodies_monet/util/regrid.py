# SPDX-License-Identifier: Apache-2.0
#

"""Model-agnostic regridder.

Single ``regrid(src, target, method=...)`` entry point that dispatches on
source/target geometry. 

First-pass implementation supports unstructured
sources sampled at arbitrary target points via a cached :class:`cKDTree`
(``nearest_s2d`` / ``nearest_d2s`` / within-radius mean).

The intent is to keep all per-model regridding logic out of model
readers and out of analysis layer call sites. Each reader just
returns an :class:`xr.Dataset` / :class:`ux.UxDataset`; 

callers that want to move values between grids call :func:`regrid` and the dispatcher
picks the right backend.

Future backends (xregrid for conservative, xesmf for structured-to-
structured, etc.) plug in without changing this API.

Usage
-----
Forward (sample unstructured model at swath pixel locations):

>>> from melodies_monet.util.regrid import regrid
>>> out = regrid(
...     modobj,
...     target={"lon": obs["lon"].values, "lat": obs["lat"].values},
...     method="nearest_s2d",
...     target_dims=("x", "y"),
... )

Backward (project flattened swath onto unstructured model columns):

>>> out = regrid(
...     flattened_swath,
...     target={"lon": mod_lon_1d, "lat": mod_lat_1d},
...     method="nearest_s2d",
...     target_dims=("n_face",),
... )
"""

import numpy as np
import xarray as xr

_LON_NAMES = ("longitude", "lon", "Longitude")
_LAT_NAMES = ("latitude", "lat", "Latitude")

def regrid(
    src,
    target,
    method="nearest_s2d",
    src_lon=None,
    src_lat=None,
    target_dims=None,
    max_distance=None,
    radius=None,
    backend="auto",
):
    """Model-agnostic regrid dispatcher.

    Parameters
    ----------
    src : xr.Dataset or ux.UxDataset
        Source data. Must carry 1-D longitude/latitude coords identifying
        its "column" dim. Other dims (time, z, ...) and attrs are
        preserved through the regrid.
    target : dict
        Target spec. Must contain ``"lon"`` and ``"lat"`` array-likes of
        the same shape (1-D, 2-D, or higher). Output preserves that
        shape in its new spatial dim(s).
    method : str, default "nearest_s2d"
        Regrid method. Currently supported:

        - ``"nearest_s2d"`` / ``"nearest_d2s"``: nearest-neighbor lookup
          from source to each target point.
        - ``"radius_mean"``: mean of all sources within ``radius`` of
          each target. Requires ``radius`` kwarg.

        Future: ``"conservative"``, ``"bilinear"`` (via xregrid).
    src_lon, src_lat : str, optional
        Names of source longitude/latitude coords. Auto-detected from
        ``longitude``/``lon``/``Longitude`` (and similarly for lat) if
        not given.
    target_dims : tuple of str, optional
        Names for the output spatial dim(s). Length must match the
        ndim of target lon/lat. Defaults: ``("target",)`` for 1-D
        target, ``("dim_0", "dim_1", ...)`` for N-D.
    max_distance : float, optional
        Nearest-neighbor distance cap (lon/lat Euclidean degrees).
        Targets with no source inside this radius NaN out.
    radius : float, optional
        Averaging radius for ``method="radius_mean"`` (degrees).
    backend : str, default "auto"
        Force a specific backend or let dispatcher choose. Currently
        only ``"ckdtree"`` is implemented.

    Returns
    -------
    xr.Dataset
        Same data_vars as ``src`` but with the source column dim
        replaced by ``target_dims``, and ``longitude``/``latitude``
        attached as coords on those dims.
    """
    src_lon = _resolve_coord_name(src, src_lon, _LON_NAMES, "longitude")
    src_lat = _resolve_coord_name(src, src_lat, _LAT_NAMES, "latitude")

    if not isinstance(target, dict) or "lon" not in target or "lat" not in target:
        raise TypeError(
            "regrid: target must be a dict with 'lon' and 'lat' array-likes. "
            f"Got: {type(target).__name__}"
        )
    tlon = np.asarray(target["lon"])
    tlat = np.asarray(target["lat"])
    if tlon.shape != tlat.shape:
        raise ValueError(
            f"regrid: target lon shape {tlon.shape} != lat shape {tlat.shape}"
        )

    if backend in ("auto", "ckdtree"):
        return _regrid_ckdtree(
            src, src_lon, src_lat, tlon, tlat, target_dims,
            method=method, max_distance=max_distance, radius=radius,
        )
    raise NotImplementedError(
        f"regrid: backend={backend!r} not implemented. Available: 'ckdtree'."
    )


def _resolve_coord_name(obj, given, candidates, kind):
    """Find a coord by name. Use ``given`` if set, else first candidate present."""
    if given is not None:
        if given not in obj.variables:
            raise KeyError(
                f"regrid: requested {kind} name {given!r} not in dataset. "
                f"Available: {list(obj.variables)}"
            )
        return given
    for c in candidates:
        if c in obj.variables:
            return c
    raise KeyError(
        f"regrid: no source {kind} found. Looked for {candidates}. "
        f"Available: {list(obj.variables)}"
    )


def _regrid_ckdtree(
    src, src_lon, src_lat, tlon, tlat, target_dims,
    method="nearest_s2d", max_distance=None, radius=None,
):
    """cKDTree backend.

    Reuses :func:`melodies_monet.util.uxarray_util.sample_unstructured_at_points`,
    which carries the cached KDTree (the same source grid hit across many
    target queries — e.g. per-granule TEMPO loops — pays the build cost
    once).
    """
    from melodies_monet.util.uxarray_util import sample_unstructured_at_points

    # The helper hardcodes "longitude"/"latitude" internally. Rename src's
    # coords if they use a different convention; restore on exit.
    rename_in = {}
    if src_lon != "longitude":
        rename_in[src_lon] = "longitude"
    if src_lat != "latitude":
        rename_in[src_lat] = "latitude"
    src_norm = src.rename(rename_in) if rename_in else src

    target_shape = tlon.shape
    tlon_flat = tlon.ravel()
    tlat_flat = tlat.ravel()

    if method == "radius_mean":
        if radius is None:
            raise ValueError(
                "regrid: method='radius_mean' requires radius=<float deg>."
            )
        sampled = sample_unstructured_at_points(
            src_norm, tlon_flat, tlat_flat, radius=radius,
        )
    elif method in ("nearest_s2d", "nearest_d2s"):
        sampled = sample_unstructured_at_points(
            src_norm, tlon_flat, tlat_flat, max_distance=max_distance,
        )
    else:
        raise NotImplementedError(
            f"regrid (ckdtree backend): method={method!r} not supported. "
            "Use 'nearest_s2d', 'nearest_d2s', or 'radius_mean'."
        )

    # Resolve target_dims.
    if target_dims is None:
        if len(target_shape) == 1:
            target_dims = ("target",)
        else:
            target_dims = tuple(f"dim_{i}" for i in range(len(target_shape)))
    target_dims = tuple(target_dims)
    if len(target_dims) != len(target_shape):
        raise ValueError(
            f"regrid: target_dims has length {len(target_dims)} but target "
            f"shape is {target_shape}."
        )

    # Fast path: 1-D target with default name -- sample_unstructured_at_points
    # already returns exactly this shape.
    if len(target_shape) == 1 and target_dims == ("target",):
        return sampled

    # 1-D target with caller-chosen dim name: just rename.
    if len(target_shape) == 1:
        return sampled.rename({"target": target_dims[0]})

    # N-D target: reshape every var's trailing "target" axis back to
    # target_shape, then attach lon/lat as coords on the target dims.
    out = xr.Dataset(attrs=dict(sampled.attrs))
    for v in sampled.data_vars:
        da = sampled[v]
        if "target" not in da.dims:
            out[v] = da
            continue
        transposed = da.transpose(..., "target")
        arr = np.asarray(transposed.values)
        new_shape = arr.shape[:-1] + target_shape
        new_dims = transposed.dims[:-1] + target_dims
        out[v] = xr.DataArray(
            arr.reshape(new_shape), dims=new_dims, attrs=dict(da.attrs),
        )

    out = out.assign_coords({
        "longitude": (target_dims, tlon),
        "latitude": (target_dims, tlat),
    })
    return out