
"""Validate the TROPOMI averaging-kernel operator on ONE granule + MPAS.

Check the AK physics and the regrid on
TROPOMI's (y, x)/longitude/bounds conventions before wiring into MM.

"""

import numpy as np
import xarray as xr

# --------------------------------------------------------------------------
TROPOMI = "/glade/campaign/acom/acom-weather/amirrezaei/tropomi_no2/2024/S5P_OFFL_L2__NO2____20240101T074458_20240101T092629_32219_03_020600_20240102T234600.nc"
MODEL_GLOB = (
    "/glade/campaign/acom/acom-weather/wenfut/MPAS_ASIAAQ/model_output/"
    "ASIA-AQ-MPAS_2024_x20.835586.grid_asiaaq_g17_mg17_58L/atm/hist/"
    "ASIA-AQ-MPAS_2024_x20.835586.grid_asiaaq_g17_mg17_58L.cam.h2i.2024-01-01-*.nc"
)
MESH = "/glade/work/wenfut/MPAS_tools/NCL_codes2/x20.835586.real.asiaaq.init_58L.nc"
MOL_M2_TO_MOLEC_CM2 = 6.02214e19
# --------------------------------------------------------------------------



def main():
    import glob
    import monetio as mio
    from melodies_monet.util.regrid_util import regrid
    from melodies_monet.util.uxarray_util import uxgrid_from_corner_bounds

    # tropo generic reader
    vd = {
        "nitrogendioxide_tropospheric_column": {},
        "averaging_kernel": {},
        "air_mass_factor_troposphere": {},
        "air_mass_factor_total": {},
        "tm5_tropopause_pressure": {},
        "latitude_bounds": {},
        "longitude_bounds": {},
        "qa_value": {},
    }
    obs = mio.sat.tropomi_l2.open_datasets(TROPOMI, vd)
    o = obs[list(obs)[0]].squeeze("time")     # drop the length-1 time dim
    print("TROPOMI dims:", dict(o.sizes))      # (z=34, y=4173, x=450, ...)
    olon = np.asarray(o["longitude"].values)   # (y, x)
    olat = np.asarray(o["latitude"].values)
    ny, nx = olon.shape
    print(f"swath (y,x) = {olon.shape}")

    # mpas model
    m = mio.models._cam_unstructured.open_mfdataset(
        sorted(glob.glob(MODEL_GLOB))[:1],
        var_list=["NO2"], mesh_file=MESH, convert_to_ppb=True,
    ).squeeze("time", drop=True)
    # Mimic _model.py: reader emits lon/lat; downstream code wants
    # longitude/latitude as coords.
    _ren = {}
    if "lon" in m.variables and "longitude" not in m.variables:
        _ren["lon"] = "longitude"
    if "lat" in m.variables and "latitude" not in m.variables:
        _ren["lat"] = "latitude"
    if _ren:
        m = m.rename(_ren)
    m = m.set_coords([c for c in ("longitude", "latitude") if c in m.variables])

    no2_units = str(m["NO2"].attrs.get("units", "")).lower()
    print("model dims:", dict(m.sizes), "| NO2 units:", no2_units or "(none)")
    # Partial-column formula expects the mixing ratio in mol/mol. NO2 is
    # either ppbV (reader converted) or mol/mol (not converted) -> pick the
    # factor accordingly so we don't get a 1e9 error.
    mr_to_molmol = 1e-9 if "ppb" in no2_units else 1.0
    print(f"  mixing-ratio->mol/mol factor: {mr_to_molmol:g}")
    col_dim = m["longitude"].dims[0]
    nlev = m.sizes["z"]

    # regrid model (NO2, pres_pa_mid, dz_m, temperature_k) to swath 
    src_vars = ["NO2", "pres_pa_mid", "dz_m", "temperature_k"]
    msrc = m[src_vars]
    target = {"lon": olon, "lat": olat}
    print("regridding model -> TROPOMI swath (nearest)...")
    on_swath = regrid(msrc, target=target, method="nearest_s2d",
                      target_dims=("y", "x"))
    # on_swath has dims (z, y, x) for each var
    print("model-on-swath dims:", dict(on_swath.sizes))

    # model partial column per MODEL layer (molec/cm2)
    no2_mod = (on_swath["NO2"] * mr_to_molmol)      # -> mol/mol
    p_mod = on_swath["pres_pa_mid"]                 # (z=58, y, x), Pa
    p_trop = o["pres_pa_mid"]                       # (z=34, y, x), Pa
    nz_t = p_trop.sizes["z"]
    no2_t = np.full((nz_t, ny, nx), np.nan, dtype=float)
    pm = np.asarray(p_mod.transpose("z", "y", "x").values)
    nm = np.asarray(no2_mod.transpose("z", "y", "x").values)
    pt = np.asarray(p_trop.transpose("z", "y", "x").values)
    for j in range(ny):
        for i in range(nx):
            mp = pm[:, j, i]; mc = nm[:, j, i]; tp = pt[:, j, i]
            good = np.isfinite(mp) & np.isfinite(mc)
            if good.sum() < 2 or not np.isfinite(tp).any():
                continue
            order = np.argsort(mp[good])
            no2_t[:, j, i] = np.interp(
                np.log10(tp), np.log10(mp[good][order]), mc[good][order],
                left=np.nan, right=np.nan,
            )
    no2_t = xr.DataArray(no2_t, dims=("z", "y", "x"))     # mol/mol on TROPOMI layers

    # partial column on the tropomi layers 
    g, M_air, NA = 9.80665, 0.0289644, 6.022e23
    pint = o["pres_pa_int"].transpose("z_stagg", "y", "x")   # (35, y, x), Pa
    dp = np.abs(
        pint.isel(z_stagg=slice(0, -1)).values
        - pint.isel(z_stagg=slice(1, None)).values
    )                                                        # (34, y, x)
    dp = xr.DataArray(dp, dims=("z", "y", "x"))
    pcol_t = no2_t * dp * (NA / (g * M_air) / 1e4)           # molec/cm2 per layer

    # AK to tropospheric AK, sum over troposphere 
    ak = o["averaging_kernel"].transpose("z", "y", "x")
    ak_trop = ak * (o["air_mass_factor_total"] / o["air_mass_factor_troposphere"])
    trop = p_trop.transpose("z", "y", "x") >= o["tm5_tropopause_pressure"]
    model_col = (ak_trop * pcol_t).where(trop).sum("z", skipna=True)   # molec/cm2
    model_col = model_col.where(np.isfinite(pcol_t.isel(z=0)))

    # compare
    obs_col = o["nitrogendioxide_tropospheric_column"] * MOL_M2_TO_MOLEC_CM2
    qa = o["qa_value"]
    keep = (qa >= 0.5)
    mv = np.asarray(model_col.where(keep).values).ravel()
    ov = np.asarray(obs_col.where(keep).values).ravel()
    ok = np.isfinite(mv) & np.isfinite(ov) & (ov > 0) & (mv > 0)
    print(f"\npaired (qa>=0.5): {ok.sum()} pixels")
    print(f"obs   median: {np.median(ov[ok]):.3e} molec/cm2")
    print(f"model median: {np.median(mv[ok]):.3e} molec/cm2")
    print(f"model/obs median ratio: {np.median(mv[ok])/np.median(ov[ok]):.2f}")
    print("  (TROPOMI trop NO2 ~1e15 clean to ~1e16 polluted; ratio ~0.5-3 = sane)")


if __name__ == "__main__":
    main()
