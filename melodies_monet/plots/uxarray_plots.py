# SPDX-License-Identifier: Apache-2.0
#
"""Unstructured-grid plot functions.

Mirror of :mod:`melodies_monet.plots.xarray_plots`: one function per plot
type, but for data on a uxarray-backed unstructured grid (e.g. CESM-SE
columns over ``n_face`` / ``ncol``). Structured-grid plots stay in
``xarray_plots``; each function there dispatches here when the input
``dset`` is unstructured.

All functions ultimately render through
:func:`melodies_monet.plots.uxarray_render.render_unstructured_field` so the
polygon construction + cartopy feature/gridline/colorbar work happens in one
place.
"""
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from melodies_monet.plots import savefig
from melodies_monet.plots.uxarray_render import render_unstructured_field


def _resolve_uxgrid(dset, uxgrid):
    """Return a uxarray ``Grid``, opening from ``dset`` attrs if needed."""
    if uxgrid is not None:
        return uxgrid
    grid_file = dset.attrs.get("mio_scrip_file") or dset.attrs.get("mio_grid_file")
    if not grid_file:
        raise ValueError(
            "uxarray plot called on unstructured data but no uxgrid was "
            "passed and dset has no mio_scrip_file/mio_grid_file attr."
        )
    import uxarray as ux

    return ux.open_grid(grid_file)


def _domain_extent(dset, domain_type, domain_name):
    """Replicate the extent-resolution logic from xarray_plots so the two
    sides produce visually identical maps on the same YAML config."""
    from monet.util.tools import get_epa_region_bounds as get_epa_bounds
    from monet.util.tools import get_giorgi_region_bounds as get_giorgi_bounds

    if domain_type == "all" and domain_name == "CONUS":
        return 25.0, -130.0, 50.0, -60.0, domain_name + ": "
    if domain_type == "epa_region" and domain_name is not None:
        latmin, lonmin, latmax, lonmax, _ = get_epa_bounds(index=None, acronym=domain_name)
        return latmin, lonmin, latmax, lonmax, "EPA Region " + domain_name + ": "
    if domain_type == "giorgi_region" and domain_name is not None:
        latmin, lonmin, latmax, lonmax, _ = get_giorgi_bounds(index=None, acronym=domain_name)
        return latmin, lonmin, latmax, lonmax, "Giorgi Region " + domain_name + ": "
    if domain_name == "model":
        latmin = float(dset["latitude"].min())
        latmax = float(dset["latitude"].max())
        lonmin = float(dset["longitude"].min())
        lonmax = float(dset["longitude"].max())
        return latmin, lonmin, latmax, lonmax, ""
    return -90, -180, 90, 180, (domain_name + ": ") if domain_name else ""


def make_spatial_bias_gridded(
    dset,
    varname_o=None,
    label_o=None,
    varname_m=None,
    label_m=None,
    ylabel=None,
    vdiff=None,
    nlevels=None,
    proj=None,
    outname="plot",
    domain_type=None,
    domain_name=None,
    fig_dict=None,
    text_dict=None,
    uxgrid=None,
    debug=False,
    **kwargs,
):
    """Spatial bias (model - obs) on an unstructured model grid.

    Drop-in unstructured counterpart to
    :func:`melodies_monet.plots.xarray_plots.make_spatial_bias_gridded` —
    same signature plus a ``uxgrid`` kwarg. Renders the diff as filled
    grid-cell polygons via uxarray instead of ``pcolormesh``.
    """
    if not debug:
        plt.ioff()

    fig_dict = fig_dict or {}
    text_kwargs = {**dict(fontsize=20), **(text_dict or {})}

    if ylabel is None:
        ylabel = varname_o
        if "units" in dset[varname_o].attrs:
            ylabel = f"{ylabel} ({dset[varname_o].attrs['units']})"

    # model - obs, time-mean if there's still a time dim
    diff = (dset[varname_m] - dset[varname_o]).squeeze()
    if "time" in diff.dims:
        diff = diff.mean("time")

    latmin, lonmin, latmax, lonmax, title_add = _domain_extent(dset, domain_type, domain_name)

    # symmetric diverging color scale
    if vdiff is None:
        vdiff = float(np.max((
            np.abs(diff.quantile(0.99)),
            np.abs(diff.quantile(0.01)),
        )))
    if nlevels is None:
        nlevels = 21
    clevel = np.linspace(-vdiff, vdiff, nlevels)
    cmap = mpl.cm.get_cmap("RdBu_r", nlevels - 1)
    norm = mpl.colors.BoundaryNorm(clevel, ncolors=cmap.N, clip=False)

    uxgrid = _resolve_uxgrid(dset, uxgrid)

    _proj = proj if proj is not None else ccrs.PlateCarree()
    fig = plt.figure(figsize=fig_dict.get("figsize", [15, 8]))
    ax = fig.add_subplot(1, 1, 1, projection=_proj)

    render_unstructured_field(
        ax, diff, uxgrid,
        cmap=cmap, norm=norm,
        extent=[lonmin, lonmax, latmin, latmax],
        states=fig_dict.get("states", True),
        cbar_label=r"$\Delta$" + ylabel,
        text_kwargs=text_kwargs,
    )

    timestamps = (
        f" {dset['time'][0].values.astype(str)[:16]}$-$"
        f"{dset['time'][-1].values.astype(str)[:16]}"
    )
    plt.title(
        title_add + label_m + " - " + label_o + timestamps,
        fontweight="bold", **text_kwargs,
    )
    plt.tight_layout(pad=0)
    savefig(
        outname + ".png", loc=4, logo_height=100,
        bbox_inches="tight", dpi=150,
    )
    return ax
