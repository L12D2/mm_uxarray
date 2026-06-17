
"""
track_pairing.py
================

Pairing utilities for **along-track model output** — CAM "1s_1pt" history
streams that sample the model at an aircraft's lat/lon/time directly

Needs to be paired via time-matching + vertical interpolation to the obs
pressure

this file is specific to cesm 

"""

import json

import numpy as np
import pandas as pd
import xarray as xr

def pair_track_model(model_file, obs_file, mapping, resample="60s",
                     obs_scale=None, fill_below=-9999.0, gap_warn_hours=1.0):
    """Pair a CAM '1s_1pt' along-track model file with aircraft obs.

    Parameters
    ----------
    model_file : str
        Path to the along-track model file (
            dims ``ncol``, ``lev``
            per-sample ``time``/``lat``/``lon`` coords and ``PMID`` in Pa).
    obs_file : str
        Path to the aircraft obs (netCDF). Needs time, pressure_obs (Pa), altitude, lat/lon, and
        the mapped obs variables.
    mapping : dict
        model_var: obs_var
    resample : str
        Pandas resample rule applied to the obs before pairing
    obs_scale : dict, optional
        Per-obs-var unit scaling applied after fill-masking
    fill_below : float
        Obs values ``<= fill_below`` are set to NaN Default -9999.
    gap_warn_hours : float
        Warn if obs file doesn't match the model's flight day

    Returns
    -------
    pandas.DataFrame
        One row per obs point
        Model gases in mol/mol are auto-scaled to ppb.
    """
    
    obs_scale = obs_scale or {"NO2": 0.001}

    # model ncol to time, (time, lev)
    m = xr.open_dataset(model_file, decode_times=True).swap_dims({"ncol": "time"}).sortby("time")
    pmid = m["PMID"].values                          # (time, lev) Pa
    mtime = pd.to_datetime(m["time"].values)

    # obs mask fills before resample/scale
    odf = xr.open_dataset(obs_file).to_dataframe().reset_index()
    odf["time"] = pd.to_datetime(odf["time"])
    for c in odf.select_dtypes("number").columns:
        odf.loc[odf[c] <= fill_below, c] = np.nan
    odf = (odf.set_index("time").resample(resample).mean(numeric_only=True)
              .dropna(subset=["pressure_obs"]).reset_index())
    for v, s in obs_scale.items():
        if v in odf:
            odf[v] = odf[v] * s

    # nearest model track sample per obs time
    mi = mtime.values.astype("datetime64[ns]").astype("int64")
    oi = odf["time"].values.astype("datetime64[ns]").astype("int64")
    pos = np.clip(np.searchsorted(mi, oi), 1, len(mi) - 1)
    j = np.where((oi - mi[pos - 1]) <= (mi[pos] - oi), pos - 1, pos)

    gap = np.abs(mi[j] - oi).max() / 1e9
    if gap > gap_warn_hours * 3600:
        print(f"WARNING: max obs-model time gap = {gap / 3600:.1f} h "
              f"-- wrong obs file for this model day?")

    out = {"time": odf["time"].values,
           "latitude": odf.get("latitude", odf.get("lat")).values,
           "longitude": odf.get("longitude", odf.get("lon")).values,
           "altitude": odf["altitude"].values,
           "pressure_obs": odf["pressure_obs"].values}

    for mod_v, obs_v in mapping.items():
        if mod_v not in m:
            continue
        da = m[mod_v]
        scale = 1e9 if "mol/mol" in da.attrs.get("units", "").lower() else 1.0  # -> ppb
        prof = da.values * scale
        mvals = np.full(len(odf), np.nan)
        for i, (jj, p) in enumerate(zip(j, odf["pressure_obs"].values)):
            order = np.argsort(pmid[jj])
            mvals[i] = np.interp(p, pmid[jj][order], prof[jj][order],
                                 left=np.nan, right=np.nan)
        out[obs_v] = odf[obs_v].values if obs_v in odf else np.nan
        out[f"{obs_v}_mod"] = mvals
    return pd.DataFrame(out)

def save_paired_mm(df, mapping, obs_label, model_label, out_nc):
    """Write track-paired DataFrame to be compatible with MM's aircraft paired-file 

    Produces dims (time, x) 
    
    read_analysis can reconstruct the pair
    """
    
    obs_vars = list(mapping.values())
    model_vars = list(mapping.keys())
    obs_names = set(obs_vars)
    n = len(df)
    dims = ("time", "x")
    f64 = lambda a: (dims, np.asarray(a, float).reshape(n, 1))
    f32 = lambda a: (dims, np.asarray(a, np.float32).reshape(n, 1))

    ds = xr.Dataset()
    ds["latitude"] = f64(df["latitude"])
    ds["longitude"] = f64(df["longitude"])
    ds["pressure_obs"] = f32(df["pressure_obs"])
    ds["altitude"] = f64(df["altitude"])

    for mod_v, obs_v in mapping.items():
        if obs_v in df:
            ds[obs_v] = f64(df[obs_v])
        mcol = f"{mod_v}_new" if mod_v in obs_names else mod_v
        ds[mcol] = f32(df[f"{obs_v}_mod"])

    t = pd.to_datetime(df["time"].values)
    ds = ds.assign_coords(time=("time", t))
    t0 = t[0]
    ds["time"].encoding = {"units": f"minutes since {t0:%Y-%m-%d %H:%M:%S}",
                           "calendar": "proleptic_gregorian", "dtype": "int64"}

    meta = {"type": "aircraft", "radius_of_influence": None,
            "obs": obs_label, "model": model_label,
            "model_vars": model_vars, "obs_vars": obs_vars,
            "filename": f"{obs_label}_{model_label}.nc"}
    ds.attrs = {"title": "", "format": "NetCDF-4",
                "dict_json": json.dumps(meta, indent=4),
                "group_name": f"{obs_label}_{model_label}"}
    ds.to_netcdf(out_nc)
    print("wrote", out_nc)

