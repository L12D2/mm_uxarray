# SPDX-License-Identifier: Apache-2.0
#
"""Shared uxarray rendering for unstructured model fields.

A single entry point, :func:`render_unstructured_field`, that every spatial
plot type (overlay, spatial_dist, gridded bias, satellite overlay, ...) can
call to draw a 1-D unstructured column field as filled grid-cell polygons on a
cartopy axis. This is the uxarray-native replacement for the legacy
``Plot_2D`` SE path.

Color scaling is controlled by the caller via ``norm``/``cmap`` (so linear,
log via ``SymLogNorm``, or diff via a diverging norm/cmap all work without
special-casing here). Data binding/alignment lives in
:mod:`melodies_monet.util.uxarray_util`.

take 1D unstructured fields; a uxarray mesh (grid file); plotting metadata to plot:

- polygon collection that represents the mesh cells 
- render directly on cartopy map axis 

"""

# https://uxarray.readthedocs.io/en/latest/user-guide/grid-formats.html 

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from melodies_monet.util.uxarray_util import uxda_from_columns

def render_unstructured_field(
    ax,
    field,
    uxgrid,
    *,
    cmap,
    norm=None,
    vmin=None,
    vmax=None,
    extent=None,
    coast=True,
    borders=True,
    states=True,
    gridlines=True,
    colorbar=True,
    cbar_label=None,
    cbar_kwargs=None,
    text_kwargs=None,
    periodic_elements="ignore",):
    
    """Draw a 1-D unstructured column ``field`` as polygons on a cartopy ``ax``.

    Parameters
    ----------
    ax : cartopy GeoAxes
        Axis (already created with a projection) to draw on.
    field : xarray.DataArray
        1-D model field over the unstructured column dimension, with 1-D
        ``longitude``/``latitude`` coordinates (e.g. a time-mean slice).
    uxgrid : uxarray.Grid
        Grid the field is bound to (model.uxgrid).
    cmap : matplotlib colormap
    norm : matplotlib.colors.Normalize, optional
        If given, controls color scaling (use ``SymLogNorm`` for log,
        ``TwoSlopeNorm``/``BoundaryNorm`` for diff). Takes precedence over
        ``vmin``/``vmax``.
    vmin, vmax : float, optional
        Used only when ``norm`` is None.
    extent : [lonmin, lonmax, latmin, latmax], optional
        Map extent in PlateCarree degrees.
    coast, borders, states, gridlines, colorbar : bool
        Toggle cartopy features / decorations.
    cbar_label : str, optional
    cbar_kwargs : dict, optional
        Overrides for ``figure.colorbar`` (default shrink/pad/extend).
    text_kwargs : dict, optional
        Passed to the colorbar label text.
    periodic_elements : {"ignore", "split", "exclude"}
        Antimeridian handling. "ignore" is correct (and fastest) for regional
        grids; use "split" for global grids.

    Returns
    -------
    matplotlib.collections.PolyCollection
        The polygon collection added to ``ax`` (e.g. to attach a colorbar to).
    """
    
    uxda = uxda_from_columns(field, uxgrid)

    poly = uxda.to_polycollection(periodic_elements=periodic_elements)
    # older uxarray returned (poly, corrected_to_gdf); newer returns poly only
    if isinstance(poly, tuple):
        poly = poly[0]

    poly.set_cmap(cmap)
    if norm is not None:
        poly.set_norm(norm)
    else:
        poly.set_clim(vmin=vmin, vmax=vmax)
    poly.set_edgecolor("face")
    poly.set_transform(ccrs.PlateCarree())
    ax.add_collection(poly)

    if coast:
        ax.coastlines(lw=0.5)
    if borders:
        ax.add_feature(cfeature.BORDERS, lw=0.5)
    if states:
        ax.add_feature(cfeature.STATES, lw=0.3)
    if extent is not None:
        ax.set_extent(extent, crs=ccrs.PlateCarree())

    if gridlines:
        gl = ax.gridlines(draw_labels=True, lw=1.0, color="black", alpha=0.5, linestyle=":")
        gl.top_labels = False
        gl.right_labels = False

    if colorbar:
        cbk = dict(shrink=0.8, pad=0.04, extend="both")
        if cbar_kwargs:
            cbk.update(cbar_kwargs)
        cbar = ax.figure.colorbar(poly, ax=ax, **cbk)
        if cbar_label is not None:
            cbar.set_label(cbar_label, fontweight="bold", **(text_kwargs or {}))

    return poly
