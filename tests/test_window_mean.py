"""Moving average of network_plot.py against the exact one, for a power that drops by ~40 decades.

    python3 tests/test_window_mean.py
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
from network_plot import window_mean

t = np.linspace(0.0, 1e13, 4000)
s, w = 1e11, 400
tt, m = window_mean(t, 1e-17 * np.exp(-t / s), w)
tau = t[w] - t[0]
exact = 1e-17 * s / tau * (np.exp(-(tt - tau) / s) - np.exp(-tt / s))
assert np.all(m > 0.0) and np.max(np.abs(m / exact - 1.0)) < 1e-3, np.max(np.abs(m / exact - 1.0))
assert np.allclose(window_mean(t, np.full_like(t, 3.0), w)[1], 3.0)
print("ok")
