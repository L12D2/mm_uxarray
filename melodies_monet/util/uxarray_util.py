# SPDX-License-Identifier: Apache-2.0
#

"""
unstructured model output binding to a uxarray grid 

these util functions only run if a grid file is provided in the YAML, preserving currently existing functionality

https://uxarray.readthedocs.io/en/latest/user-guide/grid-formats.html

Users should be able to provide EXODUS, UGRID, SCRIP, etc. This util file will handle all of those types automatically. 

UGRID; MPAS; SCRIP; EXODUS; ESMF; GEOS CS; ICON; FESOM2; HEALPix

05282026: EXODUS and SCRIP are supported; if ncol = num_element uxarray should handle automatically. If no, CKDTree is triggered for nearest neighbor interpolation. 

"""

import numpy as np
import uxarray as ux
import xarray as xr

def _coord(obj, names):
    """Return ``(name, coord)`` for the first of ``names`` present on ``obj``.

    Works for both ``xarray.Dataset`` and ``xarray.DataArray``; longitude /
    latitude may live in coords or (for a Dataset) data_vars.
    """
    for name in names:
        if name in obj.coords:
            return name, obj.coords[name]
        if name in getattr(obj, "data_vars", {}):
            return name, obj[name]
    return None, None


def _spatial_dim(obj):
    """Return ``(dim_name, lon_coord_name)`` for the unstructured column dim.

    Located via the longitude coordinate, which is 1-D along the column
    dimension for unstructured (e.g. CESM-SE ``ncol``) output.
    """
    name, lon = _coord(obj, ("longitude", "lon", "Longitude"))
    if lon is None or lon.ndim != 1:
        raise ValueError(
            "Could not locate a 1-D 'longitude' coordinate to identify the "
            "unstructured spatial dimension."
        )
    return lon.dims[0], name


def _lat_name(obj):
    name, lat = _coord(obj, ("latitude", "lat", "Latitude"))
    if lat is None:
        raise ValueError("Could not locate a 'latitude' coordinate.")
    return name


def _face_centers(uxgrid):
    """Return (lon, lat) face centers.

    Uses uxarray's built-in, vectorized face-center coordinates rather than a
    Python loop over faces (much faster on large grids). Falls back to a
    connectivity-based mean only if those properties are unavailable.
    """
    try:
        return np.asarray(uxgrid.face_lon.values), np.asarray(uxgrid.face_lat.values)
    except (AttributeError, ValueError):
        face_node = uxgrid.face_node_connectivity.values
        node_lon = uxgrid.node_lon.values
        node_lat = uxgrid.node_lat.values
        masked = np.where(face_node >= 0, face_node, 0)
        valid = (face_node >= 0).astype(float)
        counts = valid.sum(axis=1)
        center_lon = (node_lon[masked] * valid).sum(axis=1) / counts
        center_lat = (node_lat[masked] * valid).sum(axis=1) / counts
        return center_lon, center_lat


def _face_for_columns(uxgrid, mlon, mlat):
    """Nearest grid face for each model column (mirrors the Plot_2D cKDTree)."""
    from scipy.spatial import cKDTree

    center_lon, center_lat = _face_centers(uxgrid)
    mlon = np.asarray(mlon)
    mlon = np.where(mlon > 180, mlon - 360, mlon)
    flon = np.where(center_lon > 180, center_lon - 360, center_lon)

    tree = cKDTree(np.column_stack([flon, center_lat]))
    _, face_for_col = tree.query(np.column_stack([mlon, np.asarray(mlat)]))
    return face_for_col


def _bin_to_faces(obj, uxgrid, spatial_dim, face_for_col):
    """Average values on ``spatial_dim`` into faces, preserving other dims.

    Faces with no column map to NaN, matching the prior Plot_2D behavior.
    Accepts a Dataset or DataArray and returns the same type with the column
    dimension replaced by ``n_face``.
    """
    n_face = uxgrid.n_face
    return (
        obj.assign_coords(_face=(spatial_dim, face_for_col))
        .groupby("_face")
        .mean(skipna=True)
        .reindex(_face=np.arange(n_face))
        .rename({"_face": "n_face"})
    )


def model_to_uxdataset(obj, grid_file):
    """Bind a model Dataset to a uxarray Grid, aligned to the grid's faces.

    Parameters
    ----------
    obj : xarray.Dataset
        Model output on an unstructured column dimension (e.g. CESM-SE
        ``ncol``) with 1-D ``longitude``/``latitude`` coordinates.
    grid_file : str
        Path to a uxarray-readable grid file (e.g. EXODUS).

    Returns
    -------
    uxarray.UxDataset
        Dataset whose spatial dimension is the grid's ``n_face`` and whose
        ordering matches the grid faces. ``obj`` is not modified.
    """
    uxgrid = ux.open_grid(grid_file)
    spatial_dim, lon_name = _spatial_dim(obj)
    lat_name = _lat_name(obj)

    if obj.sizes[spatial_dim] == uxgrid.n_face:
        aligned = obj.rename({spatial_dim: "n_face"})
    else:
        face_for_col = _face_for_columns(
            uxgrid, obj[lon_name].values, obj[lat_name].values
        )
        spatial_vars = [v for v in obj.data_vars if spatial_dim in obj[v].dims]
        aligned = _bin_to_faces(obj[spatial_vars], uxgrid, spatial_dim, face_for_col)

    return ux.UxDataset(aligned, uxgrid=uxgrid)


def uxda_from_columns(da, uxgrid):
    """Build a face-aligned :class:`uxarray.UxDataArray` from a 1-D column field.

    ``da`` is a model field on the unstructured column dimension with 1-D
    ``longitude``/``latitude`` coordinates (e.g. a time-mean CESM-SE slice,
    possibly already region-cropped). Returns a UxDataArray over the grid's
    ``n_face`` dimension, ready for uxarray plotting. ``da`` is not modified.
    """
    spatial_dim, lon_name = _spatial_dim(da)
    lat_name = _lat_name(da)
    name = da.name if da.name is not None else "_v"
    n_face = uxgrid.n_face

    if da.sizes[spatial_dim] == n_face:
        aligned = da.rename({spatial_dim: "n_face"})
        return ux.UxDataArray(aligned, uxgrid=uxgrid)

    face_for_col = _face_for_columns(
        uxgrid, da[lon_name].values, da[lat_name].values)

    if da.ndim == 1:
        values = np.asarray(da.values, dtype=float)
        fc = np.asarray(face_for_col)
        valid = ~np.isnan(values)
        face_sum = np.bincount(fc[valid], weights=values[valid], minlength=n_face)
        face_cnt = np.bincount(fc[valid], minlength=n_face)
        binned = np.where(face_cnt > 0, face_sum / np.maximum(face_cnt, 1), np.nan)
        aligned = xr.DataArray(binned, dims=["n_face"], name=name)
        
    else:
        aligned = _bin_to_faces(
            da.to_dataset(name=name), uxgrid, spatial_dim, face_for_col
        )[name]

    return ux.UxDataArray(aligned, uxgrid=uxgrid)

# uncomment when can figure out xregrid 
# def sample_unstructured_at_points(modobj, target_lon, target_lat):
#     """Nearest-neighbor sample a model on a 1-D unstructured column dim at
#     arbitrary ``(lon, lat)`` target points.

#     Reusable beyond the satellite pipeline: any caller with a list of target
#     points (swath pixels, AirNow stations, sondes, ...) can sample an
#     unstructured model at those locations without xESMF (which OOMs on large
#     ``ncol`` sources because it treats them as degenerate structured grids).

#     Parameters
#     ----------
#     modobj : xarray.Dataset
#         Model output on a 1-D unstructured column dim with 1-D
#         ``longitude``/``latitude`` coordinates (e.g. CESM-SE after the
#         rename/promote in :mod:`melodies_monet.driver._model`).
#     target_lon, target_lat : 1-D array-like
#         Lon/lat of the target points. Either convention (0..360 or
#         -180..180) works; both sides are normalized to -180..180 before the
#         KDTree query.

#     Returns
#     -------
#     xarray.Dataset
#         Same data_vars as ``modobj`` but with the column dim replaced by a
#         1-D ``target`` dim of length ``len(target_lon)``. Other dims
#         (time, z, ...) and attrs are preserved.
#     """
    
#     from scipy.spatial import cKDTree

#     modobj = modobj.load()
    
#     mlon = np.asarray(modobj["longitude"].values)
#     mlat = np.asarray(modobj["latitude"].values)
#     mlon = ((mlon + 180.0) % 360.0) - 180.0

#     tlon = np.asarray(target_lon, dtype=float)
#     tlat = np.asarray(target_lat, dtype=float)
#     tlon = ((tlon + 180.0) % 360.0) - 180.0

#     # cKDTree.query rejects NaN/Inf inputs. Real swath data has them on edge
#     # pixels / fill values; query only finite targets and mask the rest below.
#     # valid = np.isfinite(tlon) & np.isfinite(tlat)

#     svalid = np.isfinite(mlon) & np.isfinite(mlat)
#     if not svalid.any():
#         raise ValueError(
#             "sample_unstructured_at_points: no finite source points to build "
#             "the KDTree from."
#         )
#     src_orig_idx = np.where(svalid)[0]
#     tree = cKDTree(np.column_stack([mlon[svalid], mlat[svalid]]))

#     #tree = cKDTree(np.column_stack([mlon, mlat]))
#     tvalid = np.isfinite(tlon) & np.isfinite(tlat)
#     idx = np.zeros(tlon.shape[0], dtype=np.intp)
#     if tvalid.any():
#         _, idx_in_valid = tree.query(np.column_stack([tlon[tvalid], tlat[tvalid]]))
#         idx[tvalid] = src_orig_idx[idx_in_valid]

#     col_dim = modobj["longitude"].dims[0]
#     sampled = modobj.isel({col_dim: idx}).rename({col_dim: "target"})

#     if not tvalid.all():
#         # NaN out the invalid target positions
#         sampled = sampled.where(xr.DataArray(tvalid, dims=["target"]))

#     return sampled

def sample_unstructured_at_points(
    modobj, target_lon, target_lat, max_distance=None, radius=None,):
    """Sample a model on a 1-D unstructured column dim at arbitrary
    ``(lon, lat)`` target points.

    Two aggregation modes:

    - **Nearest neighbor** (default, ``radius=None``): each target gets the
      single nearest source value. Optionally capped by ``max_distance``.
    - **Within-radius mean** (``radius=<float deg>``): each target gets the
      *mean* of all sources within ``radius`` of it. Poor-man's area-weighted
      aggregation -- use in the swath -> unstructured-model direction when
      you'd rather average all swath pixels that fall near a model cell than
      snap to the single nearest. Targets with no source in range -> NaN.

    Reusable beyond the satellite pipeline: any caller with a list of target
    points (swath pixels, AirNow stations, sondes, ...) can sample an unstructured model at those locations without xESMF.

    Parameters
    ----------
    modobj : xarray.Dataset
    Model on a 1-D unstructured column dim with 1-D ``longitude``/
    ``latitude`` coordinates.
    target_lon, target_lat : 1-D array-like
    Target points. Convention-agnostic; both sides normalized to
        ``-180..180`` before the KDTree query.
    max_distance : float, optional
        Only used when ``radius`` is None. Distance cap for nearest-neighbor;
        targets with no source within this distance NaN out.
    radius : float, optional
        If set, switches to within-radius averaging mode. Distance in
        degrees, Euclidean in lon/lat space.

    Returns
    -------
    xarray.Dataset
        Same data_vars as ``modobj`` but with the column dim replaced by a
        1-D ``target`` dim of length ``len(target_lon)``. Other dims
        (time, z, ...) and attrs are preserved.
    """

    from scipy.spatial import cKDTree

    modobj = modobj.load()
    
    mlon = np.asarray(modobj["longitude"].values)
    mlat = np.asarray(modobj["latitude"].values)
    mlon = ((mlon + 180.0) % 360.0) - 180.0

    tlon = np.asarray(target_lon, dtype=float)
    tlat = np.asarray(target_lat, dtype=float)
    tlon = ((tlon + 180.0) % 360.0) - 180.0

    # cKDTree.query rejects NaN/Inf inputs. Real swath data has them on edge
    # pixels / fill values; query only finite targets and mask the rest below.
    # valid = np.isfinite(tlon) & np.isfinite(tlat)

    svalid = np.isfinite(mlon) & np.isfinite(mlat)
    if not svalid.any():
        raise ValueError(
            "sample_unstructured_at_points: no finite source points to build "
            "the KDTree from."
        )
    src_orig_idx = np.where(svalid)[0]
    tree = cKDTree(np.column_stack([mlon[svalid], mlat[svalid]]))

    tvalid = np.isfinite(tlon) & np.isfinite(tlat)

    col_dim = modobj["longitude"].dims[0]
    n_target = tlon.shape[0]

    # === Within-radius averaging path (poor-man's area-weighted) ===
    if radius is not None:
        neighbor_lists = tree.query_ball_point(
            np.column_stack([tlon[tvalid], tlat[tvalid]]), r=radius
        )
        tvalid_pos = np.where(tvalid)[0]
        out = xr.Dataset(attrs=dict(modobj.attrs))
        for v in modobj.data_vars:
            da = modobj[v]
            if col_dim not in da.dims:
                out[v] = da
                continue
            transposed = da.transpose(..., col_dim)
            arr = np.asarray(transposed.values)
            new_dims = transposed.dims[:-1] + ("target",)
            out_arr = np.full(arr.shape[:-1] + (n_target,), np.nan, dtype=float)
            for ti, neigh in zip(tvalid_pos, neighbor_lists):
                if len(neigh) == 0:
                    continue
                sources = src_orig_idx[np.asarray(neigh, dtype=np.intp)]
                with np.errstate(invalid="ignore"):
                    out_arr[..., ti] = np.nanmean(arr[..., sources], axis=-1)
            out[v] = xr.DataArray(out_arr, dims=new_dims, attrs=dict(da.attrs))
        out = out.assign_coords(
            longitude=("target", np.asarray(tlon)),
            latitude=("target", np.asarray(tlat)),
        )
        return out

    # === Nearest-neighbor path ===
    idx = np.zeros(tlon.shape[0], dtype=np.intp)
    # Track which targets actually got a real nearest source within reach.
    # Targets without one stay False and get masked to NaN below.
    paired = np.zeros(tlon.shape[0], dtype=bool)
    
    if tvalid.any():
        query_kwargs = {}
        if max_distance is not None:
            query_kwargs["distance_upper_bound"] = max_distance
        dist, idx_in_valid = tree.query(
            np.column_stack([tlon[tvalid], tlat[tvalid]]), **query_kwargs
        )
        # cKDTree returns idx == n (i.e. len(src_orig_idx)) for "no neighbor
        # within distance_upper_bound". Treat those as unpaired.
        n_src = len(src_orig_idx)
        within = idx_in_valid < n_src
        # Avoid out-of-range indexing for unpaired entries
        safe_idx = np.where(within, idx_in_valid, 0)
        idx[tvalid] = src_orig_idx[safe_idx]
        tvalid_pos = np.where(tvalid)[0]
        paired[tvalid_pos[within]] = True
    
    sampled = modobj.isel({col_dim: idx}).rename({col_dim: "target"})
    # Mask: target keeps its sampled value only if it was both (a) finite and
    # (b) within the distance cap of an actual source point. Everything else
    # NaNs out -- this is what stops the "nearest-neighbor smear" of distant
    # swath pixels across the whole CONUS grid when only a single swath is
    # available.
    keep = tvalid & paired
    if not keep.all():
        sampled = sampled.where(xr.DataArray(keep, dims=["target"]))

    return sampled


    

