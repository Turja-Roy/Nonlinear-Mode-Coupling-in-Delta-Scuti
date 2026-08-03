import numpy as np
from scipy.integrate import cumulative_trapezoid

from coupling.background import G

INTERIOR = (0.02, 0.98)


def _interior(bg):
    return (bg.x > INTERIOR[0]) & (bg.x < INTERIOR[1])


def test_gravity(bg):
    """g = G M_r / r^2, taken via c_1 so that r = 0 is finite."""
    m = bg.r > 0
    exact = G * bg.M_r[m] / bg.r[m] ** 2
    assert np.max(np.abs(bg.g[m] - exact) / exact) < 1e-14
    assert bg.g[0] == 0.0
    # The grid ends in the atmosphere, so GM/R^2 is at x = 1, not at x[-1].
    i = int(np.argmin(np.abs(bg.x - 1.0)))
    assert np.isclose(bg.g[i], G * bg.M / bg.R**2, rtol=1e-12)


def test_enclosed_mass(bg):
    m = _interior(bg)
    num = cumulative_trapezoid(4 * np.pi * bg.r**2 * bg.rho, bg.r, initial=0.0)
    assert np.median(np.abs(num[m] - bg.M_r[m]) / bg.M_r[m]) < 1e-3


def test_dg_dr(bg):
    m = _interior(bg)
    num = np.gradient(bg.g, bg.r)
    assert np.median(np.abs(bg.dg_dr[m] - num[m]) / np.abs(num[m])) < 1e-3


def test_dlnrho_dlnr(bg):
    m = _interior(bg) & (bg.rho > 0)
    num = np.gradient(np.log(bg.rho[m]), np.log(bg.x[m]))
    assert np.median(np.abs(bg.dlnrho_dlnr[m] - num) / np.abs(num)) < 1e-3


def test_hydrostatic_equilibrium(bg):
    m = _interior(bg)
    num = np.gradient(bg.P, bg.r)
    assert np.median(np.abs(num[m] + bg.rho[m] * bg.g[m]) / np.abs(num[m])) < 1e-3


def test_dGamma1_clamped_only_in_atmosphere(bg):
    # The clamp must not reach into anything MESA actually resolved.
    assert bg.n_atm_clamped < 400
    assert bg.x[len(bg.x) - bg.n_atm_clamped] > 0.97
