"""Minimal xregrid probe: reproduce the exact-2x with two synthetic squares.

test_conservative_ones.py showed constant->2.0000 with PERFECT geometry
(pixel areas 1.000, mesh cover 1.003), on both subset and full meshes. So
the doubling is in the weight application. This strips it to the smallest
case: two unit squares -> the same two squares (identity). ones must -> 1.0.

Also dumps the Regridder's weight row-sums (expect 1.0/row; the bug -> 2.0)
and where the installed xregrid lives, plus its normalization code lines.

Tiny; safe on a login node.
"""
import inspect
import re

import numpy as np
import uxarray as ux
import xregrid
from xregrid import Regridder

print(f"xregrid: {getattr(xregrid, '__version__', '?')} at {inspect.getfile(xregrid)}",
      flush=True)

def make_grid(pad):
    """Two unit squares side by side: (0,0)-(2,1). pad>0 mimics SCRIP -1 padding."""
    node_lon = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0])
    node_lat = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    conn = np.array([[0, 1, 4, 3], [1, 2, 5, 4]], dtype=np.int64)
    if pad:
        conn = np.hstack([conn, np.full((2, pad), -1, dtype=np.int64)])
    return ux.Grid.from_topology(
        node_lon=node_lon, node_lat=node_lat,
        face_node_connectivity=conn, fill_value=-1,
    )

for pad, tag in ((0, "quads"), (6, "padded")):
    src_grid = make_grid(pad)
    tgt_grid = make_grid(0)
    src = ux.UxDataset(
        {"ones": (("n_face",), np.ones(2))}, uxgrid=src_grid)
    tgt = ux.UxDataset(
        {"_loc": (("n_face",), np.zeros(2, dtype="float32"))}, uxgrid=tgt_grid)
    rg = Regridder(src, tgt, method="conservative")
    out = rg(src)
    v = np.asarray(out["ones"].values, float)
    print(f"[{tag}] ones -> {v}  (expect [1. 1.])", flush=True)

    # weight matrix row sums, wherever xregrid stashes it
    found = False
    for attr in ("weights", "weight_matrix", "w", "A", "matrix", "_weights"):
        m = getattr(rg, attr, None)
        if m is None:
            continue
        try:
            arr = m.todense() if hasattr(m, "todense") else np.asarray(m)
            print(f"[{tag}] rg.{attr} row sums = {np.asarray(arr.sum(axis=1)).ravel()}",
                  flush=True)
            found = True
            break
        except Exception:  # noqa: BLE001
            continue
    if not found:
        print(f"[{tag}] weight attrs not found; rg attrs: "
              f"{[a for a in dir(rg) if not a.startswith('__')]}", flush=True)

# show xregrid's own normalization / apply code (line numbers for patching)
try:
    srcfile = inspect.getfile(Regridder)
    text = open(srcfile).read().splitlines()
    print(f"\n--- {srcfile}: lines mentioning norm/frac/area/apply ---", flush=True)
    for i, ln in enumerate(text, 1):
        if re.search(r"norm|frac|dstarea|def __call__|def regrid|def apply", ln, re.I):
            print(f"{i:5d}: {ln}", flush=True)
except Exception as e:  # noqa: BLE001
    print(f"source dump skipped: {e!r}", flush=True)
