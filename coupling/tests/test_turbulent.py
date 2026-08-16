"""Turbulent viscosity (Duguid+ 2020 eq. 13) and turbulent damping (MW23 eq. 11)."""

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

from coupling.modes import _viscous_F, gamma_turb, nu_turb


def _f(w):
    """The bare frequency factor, recovered from nu_turb on a one-point model."""
    class _Fake:
        v_conv = np.array([1.0])
        omega_conv = np.array([1.0])
        r = np.array([1.0])

    return float(nu_turb(_Fake(), w)[0])


def test_fit_continuous_at_breakpoints():
    for w in (1e-2, 5.0):
        assert _f(w * (1 - 1e-9)) == pytest.approx(_f(w * (1 + 1e-9)), rel=1e-6)
    assert _f(1e-3) == 5.0
    assert _f(1.0) == pytest.approx(0.5)
    assert _f(50.0) == pytest.approx(0.5 * 5**1.5 / 2500)


def test_nu_turb_lives_only_in_convection_zones(bg):
    nu = nu_turb(bg, 1e-4)
    cz = (bg.v_conv > 0) & (bg.omega_conv > 0)
    assert np.isfinite(nu).all()
    assert (nu[~cz] == 0).all()
    assert (nu[cz] > 0).all()
    # A delta Sct model has both a convective core and a thin surface zone.
    assert bg.x[cz].min() < 0.2 and bg.x[cz].max() > 0.9


def test_nu_turb_is_u_l_f(bg):
    """nu = u_mlt l_mlt f(|omega|/omega_c) on the real grid, with l_mlt =
    v_conv/omega_conv. Below the plateau it stops depending on omega."""
    cz = (bg.v_conv > 0) & (bg.omega_conv > 0)
    u, l_mlt = bg.v_conv[cz], bg.v_conv[cz] / bg.omega_conv[cz]

    omega = 2e-3  # a typical p-mode; the deep core is on the w^-2 branch here
    w = omega / bg.omega_conv[cz]
    f = np.array([_f(wi) for wi in w])
    assert nu_turb(bg, omega)[cz] == pytest.approx(u * l_mlt * f, rel=1e-12)

    lo = bg.omega_conv[cz].min() * 1e-4
    assert nu_turb(bg, lo)[cz] == pytest.approx(5.0 * u * l_mlt, rel=1e-12)
    assert nu_turb(bg, lo / 2)[cz] == pytest.approx(nu_turb(bg, lo)[cz])


def test_nu_turb_nan_without_mesa_profile(poly3):
    bg = poly3[0]
    assert np.isnan(nu_turb(bg, 1e-4)).all()


def test_viscous_F_positive_definite():
    """F is 2 eps_ij eps_ij - (2/3)(div xi)^2 after the angular integral, i.e. a
    sum of squares, so it can never go negative. This pins the bracketing of Lai
    (1994) eq. 8.8: the leading 2 multiplies dxi_r/dr alone, and the misreading
    2[(dxi_r/dr)^2 + (...)^2] gives a negative F here. It is only marginally
    positive -- the radial pieces reduce to (1/3)(2 eps_rr - eps_hh)^2, which
    vanishes on a whole family of displacements."""
    rng = np.random.default_rng(0)
    r = np.linspace(0.05, 1.0, 400)
    for l in (0, 1, 2, 3, 7):
        for _ in range(20):
            c = rng.normal(size=(2, 4))
            xi_r = sum(c[0, k] * r ** (k + 1) for k in range(4))
            xi_h = 0.0 * r if l == 0 else sum(c[1, k] * r ** (k + 1) for k in range(4))
            F = _viscous_F(l, r, xi_r, xi_h)
            assert F.min() > -1e-12 * max(abs(F).max(), 1.0)


def test_viscous_F_matches_strain_tensor_at_l0():
    """For l = 0 the angular integral is trivial and F must reduce to
    2 eps_ij eps_ij - (2/3)(div xi)^2 with eps_rr = dxi_r/dr and
    eps_tt = eps_pp = xi_r/r. This is the identity that fixes the bracketing."""
    r = np.linspace(0.05, 1.0, 400)
    for xi_r in (r**2, r**3 - 0.5 * r, np.sin(3 * r)):
        # Same derivative _viscous_F uses; the identity is exact, not numerical.
        d = CubicSpline(r, xi_r).derivative()(r)
        want = 2 * (d**2 + 2 * (xi_r / r) ** 2) - (2 / 3) * (d + 2 * xi_r / r) ** 2
        got = _viscous_F(0, r, xi_r, 0.0 * r)
        assert got == pytest.approx(want, rel=1e-12)


def test_gamma_turb_positive_and_subdominant(efs):
    """Turbulent viscosity only damps, and it stays below the radiative rate for
    the great majority of modes. MW23 quote a factor ~100; this model gives a
    median nearer 10, and a handful of weakly damped modes near the top of the
    driven band cross over, so the assertion is on the median."""
    ratio = []
    for ef in efs.values():
        g = gamma_turb(ef)
        assert np.isfinite(g) and g > 0
        if np.isfinite(ef.gamma):
            ratio.append(g / abs(ef.gamma))
    assert len(ratio) > 50
    assert np.median(ratio) < 0.5
    assert np.mean(np.array(ratio) > 1.0) < 0.1


def test_gamma_turb_scales_with_amplitude(efs):
    """gamma is quadratic in the eigenfunction, and the normalisation divides it
    out again: doubling xi at fixed E_star must quadruple the rate."""
    import dataclasses

    ef = efs[(2, 5)]
    twice = dataclasses.replace(ef, xi_r=2 * ef.xi_r, xi_h=2 * ef.xi_h,
                                div_xi=2 * ef.div_xi)
    assert gamma_turb(twice) == pytest.approx(4 * gamma_turb(ef), rel=1e-10)
