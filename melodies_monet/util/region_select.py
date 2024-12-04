"""
Masking arbitrary regions with regionmask
"""

import warnings
from ast import literal_eval

import pandas as pd

from melodies_monet.util.tools import get_epa_region_bounds, get_giorgi_region_bounds

try:
    import geopandas as gp
except ImportError:
    gp = None
try:
    import regionmask
except ImportError:
    regionmask = None

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
    except TypeError:
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


def control_custom_mask(data, domain_type, domain_name, domain_info=None, **kwargs):
    """Parses region information to return the right type of data.

    Parameters
    ----------
    data : xr.Dataset
        data to be masked
    domain_type : str
        type of data. Used to decide which function to apply. Should begin with
        "custom"
    domain_name : str
        It is used as the region name, or to read the info.
    domain_info : dict
        Dictionary containing relevant information on domain, like url, name, etc.
    **kwargs:
        Extra kwargs to pass to regionmask

    Returns
    -------
    xr.Dataset
        masked Dataset
    """
    if "custom" not in domain_type:
        raise ValueError("If regionmask is used, the domain_type should be starting with 'custom'")
    if "auto_polygon" in domain_type:
        masked_data = create_custom_mask(data, domain_info[domain_name])
    elif "defined_region" in domain_type:
        name_regiontype = domain_info[domain_name]["name_regiontype"]
        region_name = domain_info[domain_name]["region_name"]
        masked_data = create_predefined_mask(data, name_regiontype, region_name)
    elif "custom_file" in domain_type:
        mask_path = domain_info[domain_name].get("mask_path", None)
        mask_url = domain_info[domain_name].get("mask_url", None)
        region_name = domain_info[domain_name].get("region_name", None)
        masked_data = create_shapefile_mask(data, mask_path, mask_url, region_name, **kwargs)
    else:
        raise ValueError(
            "Could not identify the type of domain. Should be 'auto-polygon',"
            + " 'defined_region' or 'custom_file'"
        )
    return masked_data


def create_autoregion(data, domain_type, domain_name, domain_info=None):
    """Selects a region using predefined boundaries.

    Parameters
    ----------
    data : xr.Dataset | pd.DataFrame
        data to be masked
    domain_type : str
        type of data. Used to decide which function to apply.
        If domain_type == 'auto-region:custom_box', domain_info is required.
    domain_name : str
        This is used as the region name, or to read the info.
    domain_info: None | dict[str, tuple[float, float, float, float]]
        if not None, dict containing the domain name and a tuple with
        latmin, lonmin, latmax, lonmax. Only required if domain_type
        is auto-region:custom_box
    Returns
    -------
    xr.Dataset | pd.DataFrame
        Data as it was provided to the function
    """
    auto_region_id = domain_type.split(":")[1]
    if auto_region_id == "epa":
        bounds = get_epa_region_bounds(acronym=domain_name)
    elif auto_region_id == "giorgi":
        bounds = get_giorgi_region_bounds(acronym=domain_name)
    elif auto_region_id == "custom":
        bounds = domain_info[domain_name]
    else:
        raise ValueError(
            "Currently, region selections whithout a domain query have only "
            "been implemented for Giorgi and EPA regions. You asked for "
            f"{domain_type!r}. Soon, arbitrary rectangular boxes, US states and "
            "others will be included."
        )
    if isinstance(data, pd.DataFrame):
        data_all = data.loc[
            (data["latitude"] >= bounds[0])
            & (data["longitude"] >= bounds[1])
            & (data["latitude"] <= bounds[2])
            & (data["longitude"] <= bounds[3])
        ]
    else:
        data_all = data.where(
            (data["latitude"] >= bounds[0])
            & (data["longitude"] >= bounds[1])
            & (data["latitude"] <= bounds[2])
            & (data["longitude"] <= bounds[3]),
            drop=True,
        )
    return data_all


def select_region(data, domain_type, domain_name, domain_info=None, **kwargs):
    """Selects a region in whichever format it was provided

    Parameters
    ----------
    data : xr.Dataset | pd.DataFrame
        data to be masked
    domain_type : str
        type of data. Used to decide which function to apply.
    domain_name : str
        This is used as the region name, or to read the info.
    domain_info : dict
        Dict containing the domain_name and other relevant information, e. g.,
        lonlat box, mask_url, mask_file, etc.
    **kwargs:
        extra kwargs to pass to the selector, depending on the type of
        data or region.

    Returns
    -------
    xr.Dataset | pd.DataFrame
        Region selected. The type will be the same as the input data.
    """

    if domain_type == "all":
        return data
    if domain_type.startswith("auto-region"):
        data_masked = create_autoregion(data, domain_type, domain_name, domain_info)
    elif domain_type == "custom":
        if regionmask is None:
            raise ImportError(
                "regionmask is not installed, cannot create 'custom' type domain."
                + " If your domain is a simple box, try using auto-region:custom_box."
            )
        if domain_info is None:
            raise KeyError("If regionmask is used, domain_info must exist.")
        data_masked = control_custom_mask(data, domain_type, domain_name, domain_info, **kwargs)
    else:
        data_masked = data.query(domain_type + " == " + '"' + domain_name + '"', inplace=True)
    return data_masked
