"""
Masking arbitrary regions with regionmask
"""

import warnings
from ast import literal_eval

import geopandas as gp
import regionmask

try:
    import pooch
except ImportError:
    pooch = None


def create_custom_mask(data, mask_info):
    """Creates mask using regionmask.

    Parameters
    ----------
    data : xr.Dataset
        Dataset with the data that the user wishes to mask.
    mask_info : list[list[float, float]] | dict[str, list[list[float, float]]]
        List containing the vertices of the polygon to mask.

    Returns
    -------
    xr.Dataset
        masked data
    """
    if not isinstance(mask_info, (list, dict)):
        raise TypeError(f"mask_info type={type(mask_info)}, not valid for create_custom_mask.")
    if isinstance(mask_info, list):
        regions = regionmask.Regions([mask_info])
    else:
        all_regions = []
        abbrevs = []
        for rname, rpolygon in mask_info.items():
            abbrevs.append(rname)
            all_regions.append(rpolygon)
        regions = regionmask.Regions(all_regions, abbrevs=abbrevs, name="custom_regions")
    region_mask = regions.mask(data)
    masked_data = data.where(region_mask.notnull())
    return masked_data


def create_predefined_mask(data, name_regiontype, region_name):
    """Creates mask using regionmask.

    Parameters
    ----------
    data : xr.Dataset
        Dataset with the data that the user wishes to mask.
    name_regiontype : str
        name of regionmask defined regions (e.g., "srex", "giorgi")
    region_name : str
        region to mask. If it can be parsed to an integer, that is used.
        Otherwise, it is assumed to be an abbrev.

    Returns
    -------
    xr.Dataset
        mask for input data
    """
    regions = literal_eval(f"regionmask.defined_regions.{name_regiontype}")
    region_mask = regions.mask(data)
    try:
        selected_region = data.where(region_mask == int(region_name))
    except ValueError:
        selected_region = data.where(region_mask.cf == region_name)
    return selected_region


def create_shapefile_mask(data, mask_path=None, mask_url=None, region_name=None, **kwargs):
    """Creates mask from shapefile using regionmask and geopandas.

    Parameters
    ----------
    data : xr.Dataset
        Dataset with the data that the user wishes to mask.
    mask_path : str | pathobject
        Path to the shapefiles. They can be zipped.
    mask_url : str
        url to download the shapefiles. Requires pooch.
        If both mask_url and mask_path are provided, this will be
        ignored.
    region_name : str
        region to mask. If it can be parsed to an integer, that is used.
        Otherwise, it is assumed to be an abbrev.
    **kwargs
        Arguments to pass to regionmask.from_geopandas

    Returns
    -------
    xr.Dataset
        mask for the input data
    """

    if mask_url is not None and mask_path is not None:
        warnings.warn("mask_url and mask_path provided. Only one can be used. Selecting mask_path")

    if mask_path is not None:
        file = mask_path
    elif mask_url is not None:
        if pooch is None:
            raise ImportError("pooch is not installed, cannot download URL.")
        file = pooch.retrieve(mask_url, None)
    else:
        raise ValueError("Either mask_path or mask_url have to be provided")

    if file[-4:] == ".zip":
        file = "zip://" + file
    regions = regionmask.from_geopandas(gp.read_file(file), **kwargs)

    region_mask = regions.mask(data)
    try:
        selected_region = data.where(region_mask == int(region_name))
    except ValueError:
        selected_region = data.where(region_mask.cf == region_name)
    return selected_region
