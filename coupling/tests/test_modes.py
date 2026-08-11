"""Validation suite items 1 and 2."""

import numpy as np
import pytest

from coupling.gyre_io import load_h5
from coupling.modes import DampingRates, build_mode_list, load_model
from coupling.tests.conftest import MODEL

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


def test_gamma_lookup_reaches_the_ends():
    """The modes at either end of the list are matched on the same tolerance as
    the interior ones. Taking the mean of the two adjacent gaps halves `spacing`
    there, which is what dropped l = 6, n = +20 from the wide net."""
    w = np.arange(10.0, 20.0)  # unit spacing
    rates = DampingRates({2: w}, {2: np.arange(10.0)})
    for j in (0, 4, len(w) - 1):
        assert rates(2, w[j] + 0.2) == pytest.approx(j)  # inside tol = 0.3
        assert np.isnan(rates(2, w[j] + 0.4))


def test_gamma_lookup_ignores_a_widened_gap():
    """A gap the nonadiabatic run doubled by missing a mode must not buy extra
    tolerance -- the miss is exactly when the match is least trustworthy."""
    w = np.array([10.0, 11.0, 13.0, 14.0])  # the mode at 12 was missed
    rates = DampingRates({2: w}, {2: np.arange(4.0)})
    assert np.isnan(rates(2, 11.4))  # 0.4 > tol * 1.0, not tol * 1.5
    assert rates(2, 11.2) == pytest.approx(1.0)


def test_wide_net_gamma_coverage():
    """No false negatives on the real net. Every mode left without gamma has to
    be one the nonadiabatic run did not return -- measured against the *adiabatic*
    spacing, so the check does not just restate the matcher's own tolerance."""
    gyre = MODEL / "gyre"
    nad = ("summary_nad_wide.h5", "summary_nad_wide_hi.h5")
    if not any(gyre.glob("detail_wide/detail.*.h5")) or not all(
        (gyre / n).exists() for n in nad
    ):
        pytest.skip("wide net not built")

    _, efs = load_model(
        MODEL, detail_dir="detail_wide", inlist="gyre_ad_wide.in", nad=nad,
    )
    assert np.isfinite(efs[(6, 20)].gamma)  # 6.7e-3 nonadiabatic shift, not a miss

    w_nad: dict[int, np.ndarray] = {}
    for name in nad:
        d = load_h5(gyre / name)
        for l, w in zip(d["l"], np.real(d["omega"])):
            w_nad.setdefault(int(l), []).append(float(w))
    w_nad = {l: np.sort(np.array(v)) for l, v in w_nad.items()}

    missing = 0
    for l in sorted({l for l, _ in efs}):
        row = sorted(((ef.omega_dimless, ef) for (ll, _), ef in efs.items() if ll == l),
                     key=lambda t: t[0])
        w = np.array([wi for wi, _ in row])
        for i, (wi, ef) in enumerate(row):
            if np.isfinite(ef.gamma):
                continue
            missing += 1
            gaps = ([w[i] - w[i - 1]] if i else []) + ([w[i + 1] - w[i]] if i + 1 < len(w) else [])
            near = float(np.min(np.abs(w_nad[l] - wi)))
            assert near > 0.7 * min(gaps), (
                f"l={l} n={ef.n_pg}: nearest nad mode {near / min(gaps):.2f} spacings "
                "away -- a shift the matcher rejected, not a miss"
            )
    assert missing == 14  # 7 at the pass-1 band floor, 7 in the pass-2 g-net


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
