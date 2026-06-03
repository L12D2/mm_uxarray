
"""Standalone test: conservative regrid of an unstructured CESM-SE model
onto ONE TEMPO L2 swath, using the swath's NATIVE pixel corner bounds.

This validates `_conservative_mod2swath` (the forward model->swath
direction) before running the full `pair_data()`. It deliberately reads
TEMPO bounds straight from the L2 file (geolocation group) so it does NOT
depend on MM's obs reader carrying them 

1. Open a TEMPO L2 NO2 file: lat/lon centers + latitude_bounds/
   longitude_bounds (x, y, corner=4) from the geolocation group.
2. Build a minimal `obsobj` (lon/lat + *_bounds) and a `modobj`
   (surface NO2 on ncol, with mio_scrip_file attr).
3. Call the real `_conservative_mod2swath` and sanity-check the result.


"""

import numpy as np
import xarray as xr

# --------------------------------------------------------------------------
# EDIT THESE PATHS
# --------------------------------------------------------------------------
# A single TEMPO L2 NO2 granule (V03). Bounds live in the geolocation group.
TEMPO = "/glade/derecho/scratch/jjacdan/TEMPO/TEMPO_NO2_L2_V03_20240101T125157Z_S001G01.nc"

MODEL = (
    "/glade/campaign/acom/acom-da/conus_outputs/"
    "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.ERA5_ref_dust_M1.1.002/H1/"
    "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.ERA5_ref_dust_M1.1.002"
    ".cam.h1.2024-01-01-03600.nc"
)
SCRIP = "/glade/campaign/acom/MUSICA/grids/ne0CONUSne30x8/ne0CONUS_ne30x8_np4_SCRIP.nc"
VAR = "NO2"
METHOD = "conservative"   # or "bilinear"
# --------------------------------------------------------------------------


def load_tempo_swath(path):
    """Minimal TEMPO obsobj with lon/lat (x,y) + corner bounds (x,y,4).

    TEMPO L2 NetCDF is grouped: geolocation/ holds latitude, longitude,
    latitude_bounds, longitude_bounds. Variable/dim names can differ by
    product version -- the script prints what it finds so you can adjust.
    """
    geo = xr.open_dataset(path, group="geolocation")
    print("geolocation vars:", list(geo.variables))
    print("geolocation dims:", dict(geo.sizes))

    # Centers. TEMPO dims are typically (mirror_step, xtrack); rename to (x, y).
    lat = geo["latitude"]
    lon = geo["longitude"]
    center_dims = lat.dims
    rename = {center_dims[0]: "x", center_dims[1]: "y"}
    lat = lat.rename(rename)
    lon = lon.rename(rename)

    latb = geo["latitude_bounds"]
    lonb = geo["longitude_bounds"]
    # bounds dims: (mirror_step, xtrack, corner) -> (x, y, corner)
    b_rename = {latb.dims[0]: "x", latb.dims[1]: "y", latb.dims[2]: "corner"}
    latb = latb.rename(b_rename)
    lonb = lonb.rename(b_rename)

    obs = xr.Dataset(
        {
            "latitude_bounds": latb,
            "longitude_bounds": lonb,
        },
        coords={"lon": (("x", "y"), lon.values), "lat": (("x", "y"), lat.values)},
    )
    print(f"swath shape (x, y) = {obs['lon'].shape}, corners = {latb.sizes['corner']}")
    return obs


def load_model(path, scrip):
    """Minimal modobj: surface VAR on ncol + mio_scrip_file attr."""
    m = xr.open_dataset(path)
    da = m[VAR]
    if "time" in da.dims:
        da = da.isel(time=0)
    if "lev" in da.dims:
        da = da.isel(lev=-1)  # surface
    mod = da.to_dataset(name=VAR)
    # 1-D lon/lat coords on ncol (so regrid()'s col-dim match + coords work).
    mod = mod.assign_coords(
        longitude=("ncol", m["lon"].values),
        latitude=("ncol", m["lat"].values),
    )
    mod.attrs["mio_has_unstructured_grid"] = True
    mod.attrs["mio_scrip_file"] = scrip
    return mod


def main():
    from melodies_monet.util.regrid_util import regrid
    from melodies_monet.util.uxarray_util import uxgrid_from_corner_bounds

    print("=== loading TEMPO swath (native bounds) ===")
    obs = load_tempo_swath(TEMPO)

    print("\n=== loading model + SCRIP ===")
    mod = load_model(MODEL, SCRIP)
    print(f"model {VAR}: dims={dict(mod[VAR].sizes)}")

    # Build the swath target mesh exactly as _conservative_mod2swath does.
    olon = np.asarray(obs["lon"].values)
    nx, ny = olon.shape
    n_face_expected = nx * ny
    clon = np.asarray(obs["longitude_bounds"].values).reshape(n_face_expected, -1)
    clat = np.asarray(obs["latitude_bounds"].values).reshape(n_face_expected, -1)
    swath_grid = uxgrid_from_corner_bounds(clon, clat)
    print(f"\nswath grid: n_face={swath_grid.n_face} (expect {n_face_expected}), "
          f"n_node={swath_grid.n_node}")

    print("\n=== RAW regrid output (ground truth: where did the data land?) ===")
    raw = regrid(mod, method=METHOD, src_grid=SCRIP, target_grid=swath_grid)
    print("type:", type(raw).__name__)
    print("raw.sizes:", dict(raw.sizes))
    print("raw.coords:", list(raw.coords))
    for v in raw.data_vars:
        da = raw[v]
        print(f"  data_var {v!r}: dims={da.dims}, sizes={dict(da.sizes)}")
    # Which dim (if any) matches the FACE count (= pixels)?
    face_dim = next(
        (d for d in raw[VAR].dims if raw[VAR].sizes[d] == n_face_expected), None
    )
    node_dim = next(
        (d for d in raw[VAR].dims if raw[VAR].sizes[d] == swath_grid.n_node), None
    )
    print(f"\n{VAR} is on FACE dim? {face_dim!r}   on NODE dim? {node_dim!r}")
    if face_dim is not None:
        print("  -> GOOD: data is on faces (pixels). Reshape to (x, y) is valid.")
    elif node_dim is not None:
        print("  -> PROBLEM: xregrid placed the result on NODES, not faces. "
              "Trying the dummy-face-variable fix below.")
    else:
        print("  -> data is on neither face nor node count; inspect dims above.")

    # --- Full path through the real wrapper (now that the placeholder fix
    #     forces face output). Should return NO2 on (x, y). ---
    print("\n=== _conservative_mod2swath (real wrapper, expect (x, y)) ===")
    from melodies_monet.util.sat_l2_swath_utility_tempo import _conservative_mod2swath

    res = _conservative_mod2swath(mod, obs, method=METHOD)
    print("result dims:", dict(res.sizes))
    arr = np.asarray(res[VAR].values)
    finite = np.isfinite(arr)
    print(f"{VAR}: shape={arr.shape}, {finite.sum()}/{arr.size} finite")
    if finite.any():
        print(f"  min={np.nanmin(arr):.3e}  mean={np.nanmean(arr):.3e}  "
              f"max={np.nanmax(arr):.3e}")
    ok = tuple(res[VAR].dims[-2:]) == ("x", "y") and arr.shape[-2:] == (nx, ny)
    print(f"  -> {VAR} on (x, y) pixels? {ok}")


if __name__ == "__main__":
    main()
