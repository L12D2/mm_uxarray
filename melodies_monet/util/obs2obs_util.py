"""
Need a way to compare different obs simultaneously rather than one at a time. Because of how involved this may become, i wrote a septerate 
obs2obs util folder to keep any issues with it seperate (as much as possible) from the rest of the MM 


obs2obs: compare multiple obs AND multiple models simultaneously.

MELODIES-MONET pairing is one-obs to N-models; this util works AFTER pairing,
on the saved pair objects (an.paired), treating each (pair label, variable) as
an independent *data series*. Because the satellite obs-grid pairs share one
lat/lon grid (and model-space pairs share the model mesh)

multi-obs comparison should be pure alignment -- no new pairing needed 

Series spec (each entry in a group's ``series`` list):
    - {label: <pair label>, var: <variable>, name: <legend>, color: <opt>}
      raw values ("value" mode). Only mix series with the SAME units.
    - {label: ..., obs_var: X, mod_var: Y, name: ..., mode: bias}
      model - obs from within one pair.
    - {label: ..., obs_var: X, mod_var: Y, name: ..., mode: relative_bias}
      100*(model-obs)/obs. Dimensionless units will allow to put surface ppb,
      satellite columns, and aircraft data on one axis.

Handles pair objects that are xarray Datasets (satellite: (time,y,x),
(time,n_face), or (time,) series) and pandas DataFrames (surface/aircraft
point obs). Time matching: 'none', 'daily', or {window_local: [h0, h1]}
(solar local hour ~ UTC + lon/15; the honest way to co-locate a
geostationary full-day product with a ~13:30 LT polar orbiter).

Driver file: 

    an.read_control(); an.read_analysis()
    from melodies_monet.util import obs2obs_util
    obs2obs_util.run(an.paired, an.control_dict["obs2obs"],
                     default_outdir=an.output_dir)
                     
"""

import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# use native plotting capabilities within MM 
from melodies_monet.plots import surfplots as splots
from melodies_monet.plots import xarray_plots as xrplots

from melodies_monet.plots import savefig as _mm_savefig

def _save(outname_png):
    """Native MM save (adds the logo), same call pattern as xarray_plots."""
    try:
        _mm_savefig(outname_png, logo_height=100, bbox_inches="tight", dpi=200)
    except Exception as e:  # never lose a figure over logo trouble
        print(f"obs2obs: mm savefig failed ({e!r}); plain matplotlib save.",
              flush=True)
        plt.savefig(outname_png, dpi=200, bbox_inches="tight")

def _fnum(x):
    """PyYAML parses '5.0e16' (no '+') as a STRING, which then compares
    silently-wrong; coerce every user-supplied numeric defensively."""
    return None if x is None else float(x)

def _clip2(clip):
    """clip: [lo, hi] from YAML (float, float) or None.
       Add an additioanl way for user to clip values they may think are extraneous 
    """
    return None if clip is None else (float(clip[0]), float(clip[1]))
    
def _pair_obj(paired, label):
    if label not in paired:
        raise KeyError(
            f"obs2obs: pair label {label!r} not in paired data. "
            f"Available: {sorted(paired)}. Check read.paired.filenames.")
    return paired[label].obj

def _as_dataframe(obj):
    """Surface/aircraft pairs may carry a (Multi)Index; normalize to columns."""
    return obj.reset_index() if not isinstance(obj.index, pd.RangeIndex) else obj

def _local_hour_mask_xr(da, dset, h0, h1):
    """(time, ...) mask where solar local hour (UTC + lon/15) is in [h0, h1]."""
    if "longitude" not in dset.coords and "longitude" not in dset.variables:
        print("obs2obs: window_local requested but no longitude on this pair; "
              "skipping the window for it.", flush=True)
        return None
    hours = da["time"].dt.hour + da["time"].dt.minute / 60.0
    local = (hours + dset["longitude"] / 15.0) % 24.0
    return (local >= h0) & (local <= h1)

def _subsample(flat, max_points):
    """Deterministic strided subsample so month-scale swaths (~1e8 points)
    don't blow up the boxplot frame; quartiles are unaffected."""
    if max_points and flat.size > max_points:
        step = int(np.ceil(flat.size / max_points))
        return flat[::step]
    return flat

def _extract(spec, paired, time_match, max_points=2_000_000):
    """One series: dict(name, color, ts: pandas Series over time,
    flat: 1-D valid values, grid: time-mean (y,x) DataArray or None)."""
    label = spec["label"]
    mode = spec.get("mode", "value")
    name = spec.get("name", f"{label}:{spec.get('var', mode)}")
    obj = _pair_obj(paired, label)

    if isinstance(obj, xr.Dataset):
        if mode == "value":
            da = obj[spec["var"]]
        else:
            da_o, da_m = obj[spec["obs_var"]], obj[spec["mod_var"]]
            da = (da_m - da_o) if mode == "bias" \
                else 100.0 * (da_m - da_o) / da_o.where(da_o > 0)

        # temporary filter string
        clip = _clip2(spec.get("clip"))   # [lo, hi] physical range -> NaN outside
        if clip is not None:
            da = da.where((da >= clip[0]) & (da <= clip[1]))
            
        bounds = spec.get("bounds")
        if bounds is not None and "longitude" in obj.variables:
            w, e, s, n = bounds
            da = da.where((obj["longitude"] >= w) & (obj["longitude"] <= e)
                          & (obj["latitude"] >= s) & (obj["latitude"] <= n))

        if isinstance(time_match, dict) and "window_local" in time_match:
            h0, h1 = time_match["window_local"]
            m = _local_hour_mask_xr(da, obj, h0, h1)
            if m is not None:
                da = da.where(m)

        spatial = [d for d in da.dims if d != "time"]
        ts = da.mean(dim=spatial, skipna=True).to_series() if "time" in da.dims \
            else pd.Series(dtype=float)
        if time_match == "daily" and not ts.empty:
            ts = ts.resample("1D").mean()

        vals = da.values.ravel()
        flat = _subsample(vals[np.isfinite(vals)], max_points)

        grid = None
        if {"y", "x"} <= set(da.dims):
            grid = da.mean(dim="time", skipna=True) if "time" in da.dims else da
        return dict(name=name, color=spec.get("color"), ts=ts, flat=flat,
                    grid=grid, src=obj)

    # pandas path (surface / aircraft point pairs)
    df = _as_dataframe(obj)
    if mode == "value":
        v = df[spec["var"]]
    else:
        o, m = df[spec["obs_var"]], df[spec["mod_var"]]
        v = (m - o) if mode == "bias" else 100.0 * (m - o) / o.where(o > 0)
    
    clip = _clip2(spec.get("clip"))   # [lo, hi] physical range -> NaN outside
    if clip is not None:
        v = v.where((v >= clip[0]) & (v <= clip[1]))
        
    work = pd.DataFrame({"time": pd.to_datetime(df["time"]), "v": v})
    if "longitude" in df.columns:
        work["longitude"] = df["longitude"].values
        work["latitude"] = df["latitude"].values

    bounds = spec.get("bounds")
    if bounds is not None and "longitude" in work.columns:
        w, e, s, n = bounds
        work = work[(work.longitude >= w) & (work.longitude <= e)
                    & (work.latitude >= s) & (work.latitude <= n)]

    if isinstance(time_match, dict) and "window_local" in time_match:
        if "longitude" in work.columns:
            h0, h1 = time_match["window_local"]
            lh = (work.time.dt.hour + work.time.dt.minute / 60.0
                  + work.longitude / 15.0) % 24.0
            work = work[(lh >= h0) & (lh <= h1)]
        else:
            print(f"obs2obs: window_local skipped for {label} (no longitude).",
                  flush=True)

    ts = work.set_index("time")["v"].resample(
        "1D" if time_match == "daily" else "1h").mean().dropna()
    flat = work["v"].values
    flat = _subsample(flat[np.isfinite(flat)], max_points)
    return dict(name=name, color=spec.get("color"), ts=ts, flat=flat,
                grid=None, src=obj)


# MM plotting 

def _ts_to_dataset(ts):
    """Aligned pandas Series: 1-var Dataset xrplots.make_timeseries accepts."""
    idx = pd.DatetimeIndex(ts.index, name="time")
    da = xr.DataArray(ts.values, dims=("time",), coords={"time": idx}, name="v")
    return da.to_dataset()


def multi_timeseries(series, outname, ylabel=None, title=None, avg_window=None,
                     fig_dict=None, text_dict=None):
    """Overlay every series on one axis via xrplots.make_timeseries."""
    ax = None
    for s in series:
        if s["ts"].empty:
            print(f"obs2obs: no timeseries for {s['name']}; skipping.", flush=True)
            continue
        pdict = {"color": s["color"]} if s["color"] else {}
        ax = xrplots.make_timeseries(
            _ts_to_dataset(s["ts"]), varname="v", label=s["name"], ax=ax,
            avg_window=avg_window, ylabel=ylabel,
            plot_dict=pdict, fig_dict=fig_dict, text_dict=text_dict)
    if ax is None:
        return
    tk = {"fontsize": 16, **(text_dict or {})}
    ax.legend(fontsize=tk["fontsize"] * 0.7)
    if title:
        ax.set_title(title, fontweight="bold", fontsize=tk["fontsize"])
    plt.tight_layout()
    _save(outname + ".png")
    plt.close("all")
    print(f"obs2obs: wrote {outname}.png", flush=True)


def multi_boxplot(series, outname, ylabel=None, title=None, vmin=None,
                  vmax=None, fig_dict=None, text_dict=None, showfliers=False,
                        hline=hline,): # add outliers as an option in the yaml
    """One box per series via surfplots.make_boxplot.

    The combined frame is built with pd.concat (NOT column assignment as in
    calculate_boxplot) because our series have different lengths.
    """
    cols, label_bx = {}, []
    for s in series:
        if s["flat"].size == 0:
            print(f"obs2obs: no valid data for {s['name']}; skipping box.",
                  flush=True)
            continue
        cols[s["name"]] = pd.Series(s["flat"])
        label_bx.append(dict(color=s["color"] or "gray", linestyle="-",
                             marker="x", linewidth=1.2, markersize=6.0,
                             column=s["name"], label=s["name"]))
    if not cols:
        return
    comb_bx = pd.concat(cols, axis=1)   # union index; shorter series NaN-pad
    splots.make_boxplot(comb_bx, label_bx, ylabel=ylabel, vmin=vmin, vmax=vmax,
                        outname=outname, domain_type="all",
                        domain_name=title or "", fig_dict=fig_dict,
                        text_dict=text_dict, showfliers=showfliers, debug=False)
    print(f"obs2obs: wrote {outname} (boxplot)", flush=True)


def diff_map(series, outname, ylabel=None, title=None, vdiff=None,
             domain_name=None, fig_dict=None, text_dict=None):
    """series[0] - series[1] mean map via xrplots.make_spatial_bias_gridded."""
    a, b = series[0], series[1]
    if a["grid"] is None or b["grid"] is None:
        print("obs2obs: diff_map needs two gridded (y,x) series "
              f"({a['name']}, {b['name']}); skipping.", flush=True)
        return
    if a["grid"].shape != b["grid"].shape:
        print(f"obs2obs: diff_map grids differ {a['grid'].shape} vs "
              f"{b['grid'].shape}; regrid to a common obs grid first.",
              flush=True)
        return
    ds = xr.Dataset(
        {"s_minuend": a["grid"], "s_subtrahend": b["grid"]},
        coords={"longitude": a["src"]["longitude"],
                "latitude": a["src"]["latitude"]},
    )

    if "time" in a["src"].variables or "time" in a["src"].coords:
        _t = np.asarray(a["src"]["time"].values)
        if _t.size:
            ds = ds.assign_coords(time=("time", [_t.min(), _t.max()]))
            
    # plots s_minuend - s_subtrahend; the auto-region branch zooms to the
    # finite footprint of the difference.
    xrplots.make_spatial_bias_gridded(
        ds, varname_o="s_subtrahend", label_o=b["name"],
        varname_m="s_minuend", label_m=a["name"],
        ylabel=ylabel or f"{a['name']} - {b['name']}", vdiff=vdiff,
        outname=outname, domain_type="auto-region:box",
        domain_name=domain_name or (title or "obs2obs"),
        fig_dict=fig_dict or {"states": True, "counties": False},
        text_dict=text_dict, debug=False)
    print(f"obs2obs: wrote {outname} (diff_map)", flush=True)


# run

_PLOTTERS = {"multi_boxplot": multi_boxplot,
             "multi_timeseries": multi_timeseries,
             "diff_map": diff_map}

# spatial overlay 
# spatial bias 
# 
def run(paired, config, default_outdir="."):
    """Execute every group in the ``obs2obs:`` config against loaded pairs."""
    outdir = os.path.expandvars(config.get("output_dir", default_outdir))
    os.makedirs(outdir, exist_ok=True)
    default_tm = config.get("time_match", "none")
    max_points = int(config.get("max_points", 2_000_000))

    for gname, grp in config.get("groups", {}).items():
        dp = grp.get("data_proc") or {}
        tm = grp.get("time_match", default_tm)
        gbounds = grp.get("bounds")
        gclip = grp.get("clip")
        series = []
        for spec in grp["series"]:
            if gbounds is not None and "bounds" not in spec:
                spec = {**spec, "bounds": gbounds}
            if gclip is not None and "clip" not in spec:
                spec = {**spec, "clip": gclip}
            try:
                series.append(_extract(spec, paired, tm, max_points=max_points))
            except KeyError as e:
                print(f"obs2obs: {gname}: {e}; skipping series.", flush=True)
        if not series:
            continue
        types = grp.get("type", "multi_boxplot")
        for t in ([types] if isinstance(types, str) else types):
            fn = _PLOTTERS.get(t)
            if fn is None:
                print(f"obs2obs: unknown type {t!r} "
                      f"(use {sorted(_PLOTTERS)}).", flush=True)
                continue
            out = os.path.join(outdir, f"obs2obs_{gname}.{t}")
            kwargs = dict(ylabel=grp.get("ylabel"),
                          title=grp.get("title", gname),
                          fig_dict=grp.get("fig_kwargs"),
                          text_dict=grp.get("text_kwargs"))
            if t == "multi_boxplot":
                kwargs.update(
                    vmin=_fnum(grp.get("vmin")), vmax=_fnum(grp.get("vmax")),
                    showfliers=dp.get("showfliers",
                                      grp.get("showfliers", False)),
                    hline=_fnum(dp.get("hline", grp.get("hline"))))
            if t == "multi_timeseries":
                kwargs["avg_window"] = grp.get("avg_window")
            if t == "diff_map":
                kwargs.update(vdiff=_fnum(grp.get("vdiff")),
                              domain_name=grp.get("domain_name"))
            fn(series, out, **kwargs)

