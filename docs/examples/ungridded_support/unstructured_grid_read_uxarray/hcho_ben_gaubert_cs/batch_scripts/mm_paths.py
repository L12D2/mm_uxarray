
"""Centralized output-path scheme for the satellite pairing (mm_output_v2).

Layout:  <ROOT>/<run-slug>/<instrument>/<gridtype>/
  run-slug   descriptive MUSICA run name (see RUN_SLUG)
  instrument 'tempo' | 'tropomi'
  gridtype   'model'      unstructured model mesh
             'swath'      native L2 pixels
             'grid_0p03'  TEMPO  fixed lat/lon, 0.03 deg conservative
             'grid_0p05'  TROPOMI fixed lat/lon, 0.05 deg conservative
             'grid_0p1'   common 0.1 deg CONUS obs grid (radius_mean) + series

To add a run: add one entry to RUN_SLUG. To relocate everything: set
MM_OUTPUT_ROOT. Import `paired_dir` / `RUN_SLUG` from the pairing scripts so
the scheme lives in exactly one place.
"""

import os

ROOT = os.environ.get("MM_OUTPUT_ROOT", "/glade/work/lcthompson/mm_output_v2")

# RUN env key -> descriptive path slug
RUN_SLUG = {
    "nonbiog": "musicav0_conus_nonbiog",
    "biog":    "musicav0_conus_ref",
    "grapes":  "musicav0_conus_gra2pesv2",
    "mxcat":   "musicav0_hamaq_ref",
}


def run_slug(run):
    """Descriptive slug for a RUN key (falls back to the key itself)."""
    return RUN_SLUG.get(run, run)


def instrument_of(product):
    """'tempo' or 'tropomi' from an obs product name like 'tempo_l2_hcho'."""
    p = product.lower()
    if p.startswith("tempo"):
        return "tempo"
    if p.startswith("tropomi"):
        return "tropomi"
    raise ValueError(f"cannot infer instrument from product {product!r}")


def gridtype_of(target, restag=None):
    """Map (regrid target, resolution tag) -> gridtype directory.

    target 'model'/'swath' map directly; 'obs' maps to grid_<restag>
    (e.g. restag '0p03' -> 'grid_0p03', '0p1' -> 'grid_0p1').
    """
    t = str(target).lower()
    if t == "model":
        return "model"
    if t == "swath":
        return "swath"
    if t in ("obs"):
        if not restag:
            raise ValueError("obs gridtype needs a resolution tag")
        return f"grid_{restag}"
    raise ValueError(f"unknown regrid target {target!r}")


def paired_dir(run, instrument, gridtype, make=True):
    """<ROOT>/<run-slug>/<instrument>/<gridtype>; created unless make=False."""
    d = os.path.join(ROOT, run_slug(run), instrument, gridtype)
    if make:
        os.makedirs(d, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# QA-screen sensitivity toggles (base flag always on; cloud + SZA switchable).
# Uses MM's obs-level data_proc.filter_dict (applied in open_obs, pre-regrid).
# ---------------------------------------------------------------------------
def _envflag(name):
    return os.environ.get(name, "").strip().lower() in ("1", "on", "true", "yes")


def filter_spec(product, cloud_on, sza_on):
    """(filter_dict, screen_vars) for a product given cloud/SZA toggles.

    Base quality screen is always applied (TEMPO main_data_quality_flag == 0;
    TROPOMI qa_value >= 0.75 NO2 / 0.5 HCHO,CO). cloud/SZA add pixel screens.
    screen_vars must be present in the obs `variables:` so the reader loads them.
    """
    inst = instrument_of(product)
    fd, svars = {}, set()
    if inst == "tempo":
        fd["main_data_quality_flag"] = {"oper": "<=", "value": 0}
        svars.add("main_data_quality_flag")
        if cloud_on:
            fd["eff_cloud_fraction"] = {"oper": "<=", "value": 0.2}
            svars.add("eff_cloud_fraction")
        if sza_on:
            fd["solar_zenith_angle"] = {"oper": "<=", "value": 70}
            svars.add("solar_zenith_angle")
    else:  # tropomi
        qmin = 0.75 if product.lower().endswith("no2") else 0.5
        fd["qa_value"] = {"oper": ">=", "value": qmin}
        svars.add("qa_value")
        if sza_on:
            fd["solar_zenith_angle"] = {"oper": "<=", "value": 70}
            svars.add("solar_zenith_angle")
        if cloud_on:
            # NO2 uses cloud_fraction_crb; HCHO/CO cloud var names differ, so this
            # is best-effort for NO2 and may need per-product tuning.
            fd["cloud_fraction_crb"] = {"oper": "<=", "value": 0.2}
            svars.add("cloud_fraction_crb")
    return fd, svars


def filter_tag(cloud_on, sza_on):
    """Compact filename tag, e.g. cld1sza0 (cloud on, SZA off)."""
    return f"cld{int(bool(cloud_on))}sza{int(bool(sza_on))}"


def apply_filters(obs_block, product):
    """Inject QA screening into an obs control block from CLOUD_FILTER / SZA_FILTER
    env flags; returns the filename tag. Ensures screen vars are read and sets
    data_proc.filter_dict (overriding whatever the base control had)."""
    cloud_on, sza_on = _envflag("CLOUD_FILTER"), _envflag("SZA_FILTER")
    fd, svars = filter_spec(product, cloud_on, sza_on)
    v = obs_block.setdefault("variables", {})
    for sv in svars:
        v.setdefault(sv, {})
    obs_block.setdefault("data_proc", {})["filter_dict"] = fd
    return filter_tag(cloud_on, sza_on)