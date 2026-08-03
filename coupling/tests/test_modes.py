"""Validation suite items 1 and 2."""

import numpy as np
import pytest

from coupling.modes import build_mode_list

INTERIOR = (0.02, 0.95)


def _interior(bg):
    return (bg.x > INTERIOR[0]) & (bg.x < INTERIOR[1])


def test_normalisation(efs, bg):
    """Test 1: 2 omega^2 int rho |xi|^2 d3x == E_star."""
    err = [
        abs(2 * ef.omega**2 * ef.inertia_integral() / bg.E_star - 1.0)
        for ef in efs.values()
    ]
    assert max(err) < 1e-12, f"worst {max(err):.3e}"


def test_divergence_matches_numerical(efs, bg):
    """Test 2: div.xi = r^-2 d(r^2 xi_r)/dr - Lambda^2 xi_h / r, independently."""
    m = _interior(bg)
    worst = 0.0
    r = np.where(bg.r > 0, bg.r, 1.0)
    for key, ef in efs.items():
        num = np.gradient(bg.r**2 * ef.xi_r, bg.r) / r**2 - ef.Lambda2 * ef.xi_h / r
        scale = np.percentile(np.abs(ef.div_xi[m]), 90)
        err = float(np.median(np.abs(ef.div_xi[m] - num[m])) / scale)
        worst = max(worst, err)
        assert err < 5e-2, f"{key}: {err:.3e}"
    assert worst < 5e-2


def test_euler_route_agrees(efs, bg):
    """Both routes to div.xi must agree wherever the Euler one is defined."""
    from coupling.gyre_io import load_h5
    from coupling.modes import div_xi_from_euler, div_xi_from_lag_rho

    detail = next(iter(sorted((bg_dir(bg) / "detail").glob("detail.l[123].*.h5"))))
    d = load_h5(detail)
    if "lag_rho" not in d:
        pytest.skip("lag_rho not in detail_item_list")
    omega = float(np.real(d["omega"])) * bg.omega_dyn
    a, b = div_xi_from_lag_rho(d), div_xi_from_euler(d, bg, omega)
    m = _interior(bg)
    assert np.median(np.abs(a[m] - b[m])) / np.percentile(np.abs(a[m]), 90) < 1e-6


def test_pressure_perturbation(efs, bg):
    """Horizontal Euler equation: P' = rho omega^2 r xi_h - rho dPhi."""
    from coupling.gyre_io import load_h5

    detail = next(iter(sorted((bg_dir(bg) / "detail").glob("detail.l[123].*.h5"))))
    d = load_h5(detail)
    if "eul_P" not in d:
        pytest.skip("eul_P not in detail_item_list")
    from coupling.background import G

    omega = float(np.real(d["omega"])) * bg.omega_dyn
    xi_h = np.real(d["xi_h"]) * bg.R
    dPhi = np.real(d["eul_Phi"]) * G * bg.M / bg.R
    lhs = np.real(d["eul_P"]) * bg.P
    rhs = bg.rho * (omega**2 * bg.r * xi_h - dPhi)
    m = _interior(bg)
    assert np.median(np.abs(lhs[m] - rhs[m])) / np.percentile(np.abs(lhs[m]), 90) < 1e-6


def test_gammas_attached(efs):
    assert all(np.isfinite(ef.gamma) for ef in efs.values())
    driven = [ef for ef in efs.values() if ef.gamma < 0]
    assert len(driven) == 43
    assert sum(1 for ef in driven if ef.n_pg > 0) > len(driven) / 2


def test_mode_list(efs):
    modes = build_mode_list(efs)
    assert len(modes) == sum(2 * (2 * l + 1) for l, _ in efs)
    for mode in modes:
        assert mode.omega == mode.s * mode.ef.omega
    plus = {mo.key for mo in modes if mo.s == +1}
    assert all((n, l, m, +1) in plus for n, l, m, s in (mo.key for mo in modes))


def bg_dir(bg):
    import pathlib

    return pathlib.Path(__file__).resolve().parents[2] / "models" / "dsct_M2.0" / "gyre"
