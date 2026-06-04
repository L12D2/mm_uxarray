
"""Standalone test: can we do TRUE conservative regridding of ne0CONUS
CESM-SE output by depadding the SCRIP dual cells into a clean UGRID?

ne0CONUS SCRIP stores each GLL dual cell (a 4/6/8/10-gon) padded to 10
corners by repeating corners. The repeats are degenerate and ESMF refuses
conservative weights. Here we:

  1. Read SCRIP corners, dedup repeated corners per face -> true polygons.
  2. Build a clean uxarray Grid (variable face_node_connectivity).
  3. xregrid conservative: model (n_face) -> 0.1deg global grid.
  4. Mass-conservation check: area-weighted integral before vs after.

"""

import numpy as np
import xarray as xr
import uxarray as ux

# --------------------------------------------------------------------------
# EDIT THESE PATHS
# --------------------------------------------------------------------------
SCRIP = "/glade/campaign/acom/MUSICA/grids/ne0CONUSne30x8/ne0CONUS_ne30x8_np4_SCRIP.nc"
MODEL = (
    "/glade/campaign/acom/acom-da/conus_outputs/"
    "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.ERA5_ref_dust_M1.1.002/H1/"
    "f.e22.FCnudged.ne0CONUSne30x8_ne0CONUSne30x8_mt12.ERA5_ref_dust_M1.1.002"
    ".cam.h1.2024-01-01-03600.nc"
)
VAR = "NO2"          # any per-column variable in the model file
TARGET_RES = 0.1     # degrees for the global lat/lon target
FILL = -1            # fill value for ragged face_node_connectivity
# --------------------------------------------------------------------------


def build_clean_uxgrid_from_scrip(scrip_path):
    """Depad SCRIP dual cells and return a clean uxarray Grid."""
    s = xr.open_dataset(scrip_path)
    clat = np.asarray(s["grid_corner_lat"].values)  # (n_face, 10)
    clon = np.asarray(s["grid_corner_lon"].values)
    units = str(s["grid_corner_lat"].attrs.get("units", "")).lower()
    if "rad" in units:
        clat = np.rad2deg(clat)
        clon = np.rad2deg(clon)

    nf, mc = clat.shape
    print(f"SCRIP: {nf} faces, padded to {mc} corners, units={units!r}")

    # Global node dedup via a rounded (lon, lat) key.
    key = np.round(np.stack([clon.ravel(), clat.ravel()], axis=1), 6)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    node_lon = uniq[:, 0].astype(float)
    node_lat = uniq[:, 1].astype(float)
    inv = inv.reshape(nf, mc)
    print(f"unique nodes: {len(node_lon)}")

    # Per-face consecutive dedup (SCRIP pads by repeating corners), plus
    # drop a closing duplicate if the last kept == first.
    face_conn = np.full((nf, mc), FILL, dtype=np.int64)
    for i in range(nf):
        row = inv[i]
        keep = [int(row[0])]
        for j in range(1, mc):
            if int(row[j]) != keep[-1]:
                keep.append(int(row[j]))
        if len(keep) > 1 and keep[-1] == keep[0]:
            keep.pop()
        face_conn[i, : len(keep)] = keep

    sides = (face_conn != FILL).sum(axis=1)
    u, c = np.unique(sides, return_counts=True)
    print(f"clean side-count distribution: {dict(zip(u.tolist(), c.tolist()))}")

    grid = ux.Grid.from_topology(
        node_lon=node_lon,
        node_lat=node_lat,
        face_node_connectivity=face_conn,
        fill_value=FILL,
    )
    return grid


def main():
    from xregrid import Regridder, create_global_grid

    print("=== building clean uxgrid from SCRIP ===")
    grid = build_clean_uxgrid_from_scrip(SCRIP)
    print(grid)

    print("\n=== loading model field ===")
    m = xr.open_dataset(MODEL)
    da = m[VAR]
    # surface level, first time -> 1-D (ncol,)
    if "time" in da.dims:
        da = da.isel(time=0)
    if "lev" in da.dims:
        da = da.isel(lev=-1)  # surface (assumes lev ascending top->surface)
    da = da.rename({"ncol": "n_face"})
    assert da.sizes["n_face"] == grid.n_face, (
        f"model n_face {da.sizes['n_face']} != grid n_face {grid.n_face}"
    )
    src = ux.UxDataArray(da, uxgrid=grid).to_dataset(name=VAR)

    print("\n=== building conservative regridder (the real test) ===")
    tgt = create_global_grid(TARGET_RES, TARGET_RES)
    rg = Regridder(src, tgt, method="conservative", periodic=True)
    print("WEIGHTS BUILT OK -- conservative works on the depadded mesh.")

    out = rg(src)

    print("\n=== mass-conservation check (area-weighted integral) ===")
    # Source integral: sum(value * face_area). uxarray face_areas are on the
    # unit sphere; absolute units cancel since we compare src vs tgt the same way.
    src_area = np.asarray(grid.face_areas.values)
    src_val = np.asarray(da.values, dtype=float)
    finite = np.isfinite(src_val)
    src_int = np.nansum(src_val[finite] * src_area[finite])

    # Target integral: value * cos(lat) * dlat * dlon (proportional to cell area).
    tlat = np.deg2rad(np.asarray(tgt["lat"].values))
    out_val = np.asarray(out[VAR].values, dtype=float)
    w = np.cos(tlat)[:, None] * np.ones((1, out_val.shape[-1]))
    tgt_int = np.nansum(out_val * w)

    print(f"src area-integral  = {src_int:.6e}")
    print(f"tgt area-integral  = {tgt_int:.6e} (proportional; check the RATIO)")
    print("NOTE: src and tgt area units differ; the meaningful check is that")
    print("conservative weights row-sum to 1 (xregrid enforces) and that the")
    print("global mean is preserved. Compare src/tgt MEAN below:")
    print(f"src mean = {np.nanmean(src_val):.6e}")
    print(f"tgt mean = {np.nanmean(out_val):.6e}")


if __name__ == "__main__":
    main()
