import xarray as xr
import regionmask

def create_custom_mask(data, mask_info) -> xr.Dataset:
    """Creates a 3D mask using regionmask.

    Parameters
    ----------
    data : xr.Dataset
        Dataset with the data that the user wishes to mask.
    mask_info : list[list[float, float]] | dict[str, list[list[float, float]]]
        List containing the vertices of the polygon to mask.

    Returns
    -------
    xr.Dataset
        mask for input data
    """
    if not (isinstance(mask_info, list) or isinstance(mask_info, dict)):
        raise f"mask_info type={type(mask_info)}, not valid for create_custom_mask."
    elif isinstance(mask_info, list):
        regions = regionmask.Regions([mask_info])
    else:
        all_regions = []
        abbrevs = []
        for rname, rpolygon in mask_info.items():
            abbrevs.append(rname)
            all_regions.append(rpolygon)
        regions = regionmask.Regions(all_regions, abbrevs=abbrevs, name='custom_regions')
    return regions



