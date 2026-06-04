
"""Single-granule conservative pairing test: seams vs moire.

Runs conservative model->swath (forward) and swath->model (backward) on ONE
TEMPO granule -- no granule concatenation. If diagonal banding still appears,
it's resolution-matched moire (TEMPO ~3km vs ne0CONUS ~3km), not granule
seams. If it's clean here but striped in the full run, it's seams.

Saves two PNGs:
  single_granule_mod_on_swath.png   (forward: model regridded to swath pixels)
  single_granule_mod_on_model.png   (backward: swath values back on model cols)

Run on Casper in the env with xregrid + esmpy + uxarray + MM.
"""

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# EDIT THESE PATHS  (TEMPO = a SINGLE granule file)
# --------------------------------------------------------------------------
TEMPO = "/glade/derecho/scratch/jjacdan/TEMPO/TEMPO_NO2_L2_V03_20240101T125157Z_S001G01.nc"  # one granule
MODEL = (
    "/glade/campaign/acom/acom-da/conus_outputs/"
    "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.ERA5_ref_dust_M1.1.002/H1/"
    "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.ERA5_ref_dust_M1.1.002"
    ".cam.h1.2024-01-01-03600.nc"
)
SCRIP = "/glade/campaign/acom/MUSICA/grids/ne0CONUSne30x8/ne0CONUS_ne30x8_np4_SCRIP.nc"
VAR = "NO2"
METHOD = "conservative"
# --------------------------------------------------------------------------

def load_tempo_swath(path):
    """obsobj with lon/lat (x,y) + latitude_bounds/longitude_bounds (x,y,4)."""
    geo = xr.open_dataset(path, group="geolocation")
    lat, lon = geo["latitude"], geo["longitude"]
    cdims = lat.dims
    lat = lat.rename({cdims[0]: "x", cdims[1]: "y"})
    lon = lon.rename({cdims[0]: "x", cdims[1]: "y"})
    latb, lonb = geo["latitude_bounds"], geo["longitude_bounds"]
    bd = latb.dims
    br = {bd[0]: "x", bd[1]: "y", bd[2]: "corner"}
    latb = latb.rename(br)
    lonb = lonb.rename(br)
    obs = xr.Dataset(
        {"latitude_bounds": latb, "longitude_bounds": lonb},
        coords={"lon": (("x", "y"), lon.values), "lat": (("x", "y"), lat.values)},
    )
    return obs


def load_model(path, scrip):
    m = xr.open_dataset(path)
    da = m[VAR]
    if "time" in da.dims:
        da = da.isel(time=0)
    if "lev" in da.dims:
        da = da.isel(lev=-1)
    mod = da.to_dataset(name=VAR)
    mod = mod.assign_coords(
        longitude=("ncol", m["lon"].values),
        latitude=("ncol", m["lat"].values),
    )
    mod.attrs["mio_has_unstructured_grid"] = True
    mod.attrs["mio_scrip_file"] = scrip
    return mod


def main():
    from melodies_monet.util.sat_l2_swath_utility_tempo import (
        _conservative_mod2swath,
        _conservative_swath2mod,
    )

    obs = load_tempo_swath(TEMPO)
    mod = load_model(MODEL, SCRIP)
    olon = np.asarray(obs["lon"].values)
    olat = np.asarray(obs["lat"].values)
    print(f"granule swath (x,y) = {olon.shape}")

    # ---- FORWARD: model -> swath ----
    fwd = _conservative_mod2swath(mod, obs, method=METHOD)
    mv = np.asarray(fwd[VAR].values)
    print(f"forward model-on-swath: {np.isfinite(mv).sum()}/{mv.size} finite, "
          f"mean={np.nanmean(mv):.3e}")

    # Scatter only finite pixels (swath lon/lat have fill values at edges,
    # which pcolormesh rejects).
    fok = np.isfinite(olon) & np.isfinite(olat) & np.isfinite(mv)
    plt.figure(figsize=(9, 6))
    sc = plt.scatter(olon[fok], olat[fok], c=mv[fok], s=3, marker="s")
    plt.colorbar(sc, label=f"{VAR} (model on swath)")
    plt.title("FORWARD conservative: model -> swath (1 granule)\n"
              "fine diagonal banding here = moire (no seams in 1 granule)")
    plt.xlabel("lon"); plt.ylabel("lat")
    plt.xlim(olon[fok].min() - 1, olon[fok].max() + 1)
    plt.ylim(olat[fok].min() - 1, olat[fok].max() + 1)
    plt.savefig("single_granule_mod_on_swath.png", dpi=120, bbox_inches="tight")
    print("wrote single_granule_mod_on_swath.png")

    # ---- BACKWARD: swath -> model ----
    # Build a 'paired swath' that carries the model values AND the bounds,
    # exactly what _conservative_swath2mod expects.
    paired = xr.Dataset(
        {
            VAR: (("x", "y"), mv),
            "latitude_bounds": obs["latitude_bounds"],
            "longitude_bounds": obs["longitude_bounds"],
        },
        coords={"longitude": (("x", "y"), olon), "latitude": (("x", "y"), olat)},
    )
    back = _conservative_swath2mod(paired, mod, method=METHOD)
    col_dim = mod["longitude"].dims[0]
    bv = np.asarray(back[VAR].values)
    blon = np.asarray(back["longitude"].values)
    blat = np.asarray(back["latitude"].values)
    # Model lon is 0..360; convert to -180..180 to match the swath.
    blon = ((blon + 180.0) % 360.0) - 180.0

    # Swath bounding box (from the obs centers we paired against).
    sw_ok = np.isfinite(olon) & np.isfinite(olat)
    lon0, lon1 = olon[sw_ok].min(), olon[sw_ok].max()
    lat0, lat1 = olat[sw_ok].min(), olat[sw_ok].max()
    print(f"swath bbox: lon [{lon0:.1f}, {lon1:.1f}], lat [{lat0:.1f}, {lat1:.1f}]")

    # Keep only model columns INSIDE the swath bbox (the rest are the global
    # mesh coming back as ~0 / uncovered). This zooms to where pairing happened.
    inbox = (blon >= lon0) & (blon <= lon1) & (blat >= lat0) & (blat <= lat1)
    finite = np.isfinite(bv)
    print(f"backward: {finite.sum()}/{bv.size} finite columns globally; "
          f"{(inbox & finite).sum()} inside swath bbox")
    # How many in-box columns are ~0 (uncovered) vs real?
    inbox_vals = bv[inbox & finite]
    nz = inbox_vals[inbox_vals > 0]
    print(f"  in-box: {inbox_vals.size} cols, {nz.size} with value>0, "
          f"max={inbox_vals.max():.3e}")

    sel = inbox & finite
    # Color scale from in-box non-zero percentiles so structure is visible.
    if nz.size:
        vmin, vmax = np.percentile(nz, [2, 98])
    else:
        vmin, vmax = None, None
    plt.figure(figsize=(9, 6))
    sc = plt.scatter(blon[sel], blat[sel], c=bv[sel], s=6, marker="s",
                     vmin=vmin, vmax=vmax)
    plt.colorbar(sc, label=f"{VAR} (swath back on model cols)")
    plt.title("BACKWARD conservative: swath -> model columns (1 granule, zoomed)\n"
              "diagonal banding here = moire; clean here but striped in full = seams")
    plt.xlabel("lon"); plt.ylabel("lat")
    plt.xlim(lon0 - 1, lon1 + 1)
    plt.ylim(lat0 - 1, lat1 + 1)
    plt.savefig("single_granule_mod_on_model.png", dpi=120, bbox_inches="tight")
    print("wrote single_granule_mod_on_model.png")


if __name__ == "__main__":
    main()