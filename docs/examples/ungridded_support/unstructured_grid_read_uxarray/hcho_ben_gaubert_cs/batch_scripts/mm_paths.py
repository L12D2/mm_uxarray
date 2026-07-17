
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

_TROPOMI_CLOUD_VAR = {
    "tropomi_l2_no2":  "cloud_fraction_crb",
    "tropomi_l2_hcho": "cloud_fraction_crb",
    "tropomi_l2_co":   None,
}

_SAVE_DIAG = {
    "tempo":   ["solar_zenith_angle", "eff_cloud_fraction", "main_data_quality_flag"],
    "tropomi": ["solar_zenith_angle", "qa_value"],
}

_AMF_VAR = {
    "tempo_l2_no2":    ["amf_troposphere", "amf_total"],
    "tempo_l2_hcho":   ["amf"],
    "tropomi_l2_no2":  ["air_mass_factor_troposphere", "air_mass_factor_total"],
    "tropomi_l2_hcho": ["formaldehyde_tropospheric_air_mass_factor"],
    "tropomi_l2_co":   [],
}

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
            cvar = _TROPOMI_CLOUD_VAR.get(product)
            if cvar:
                fd[cvar] = {"oper": "<=", "value": 0.2}
                svars.add(cvar)
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
    if _envflag("NO_QA"):
        fd, svars = {}, set()
    else:
        fd, svars = filter_spec(product, cloud_on, sza_on)
    v = obs_block.setdefault("variables", {})
    for sv in svars:
        v.setdefault(sv, {})
    obs_block.setdefault("data_proc", {})["filter_dict"] = fd
    return filter_tag(cloud_on, sza_on)

def save_spec(product):
    """Obs vars to mark ``save: True`` for a product: the union of the SAVE_DIAG
    per-instrument diagnostic set (sza/cloud/qa, when SAVE_DIAG=on) and any
    explicit SAVE_VARS (comma list). SAVE_VARS is ADDITIVE, so you can layer
    e.g. AMF/AK on top of the diagnostics. Empty when neither is set.

    Note on cost: AMF vars are 2-D per-pixel (cheap, any target); the averaging
    kernel is 3-D layer-resolved """
    vs = []
    if _envflag("SAVE_DIAG"):
        inst = instrument_of(product)
        vs = list(_SAVE_DIAG.get(inst, []))
        if inst == "tropomi":
            cvar = _TROPOMI_CLOUD_VAR.get(product)
            if cvar and cvar not in vs:
                vs.append(cvar)
    for v in (s.strip() for s in os.environ.get("SAVE_VARS", "").split(",")):
        if v and v not in vs:
            vs.append(v)
    return vs


def apply_save(obs_block, product):
    """Mark diagnostic obs vars ``save: True`` (and ensure they're read) so the
    pairing carries them into the paired output. Returns True if any were set."""
    vs = save_spec(product)
    v = obs_block.setdefault("variables", {})
    for sv in vs:
        cur = v.get(sv)
        v[sv] = {**cur, "save": True} if isinstance(cur, dict) else {"save": True}
    return bool(vs)


def save_suffix():
    """Filename-tag suffix ('sv') distinguishing save-diagnostics products from
    plain filtered ones"""
    return "sv" if save_spec("tempo_l2_no2") or save_spec("tropomi_l2_no2") else ""
    