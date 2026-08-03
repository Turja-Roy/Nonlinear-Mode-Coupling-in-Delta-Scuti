"""Reading GYRE HDF5 output."""

from __future__ import annotations

import pathlib

import h5py
import numpy as np


def _decompound(v):
    v = np.asarray(v)
    if v.dtype.names == ("re", "im"):
        return v["re"] + 1j * v["im"]
    return v


def load_h5(path: pathlib.Path) -> dict[str, np.ndarray]:
    """Attributes and datasets in one flat dict, complex dtypes unpacked."""
    out: dict[str, np.ndarray] = {}
    with h5py.File(path, "r") as f:
        for k, v in f.attrs.items():
            out[k] = _decompound(v)
        for k in f:
            out[k] = _decompound(f[k][...])
    return out
