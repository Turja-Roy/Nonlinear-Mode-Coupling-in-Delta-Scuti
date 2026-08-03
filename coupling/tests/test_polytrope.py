"""M0/M1: the whole pipeline on polytropes, with no MESA input at all.

The delta Sct tests can only check kappa against itself. These add structures
whose answers are known independently: n = 3 has constant Gamma_1, so
dGamma1_dlnrho_s vanishes and the A55 bracket is exact; n = 0 has constant
density and an analytic f-mode, which pins several groups to zero.
"""

import numpy as np
import pytest

from coupling.background import G
from coupling.kappa import GROUPS, grid_convergence, kappa_abc

P3 = [(2, 5), (2, 15), (2, 6)]
G3 = [(1, -18), (1, -19), (2, -20)]
MIXED_L3 = [(1, -4), (2, 7), (3, 5)]
CASES3 = [P3, G3, MIXED_L3, [(0, 9), (0, 1), (0, 6)]]


def _efs(efs, keys):
    return [efs[tuple(k)] for k in keys]


# ---------------------------------------------------------------- n = 3


def test_no_mesa_input_needed(poly3):
    bg, efs = poly3
    assert np.ptp(bg.Gamma_1) < 1e-9
    assert np.all(bg.dGamma1_dlnrho_s == 0.0)
    assert not bg.has_thermal_structure
    with pytest.raises(ValueError):
        _ = bg.t_thermal


def test_structure(poly3):
    bg, _ = poly3
    m = bg.r > 0
    assert np.max(np.abs(bg.g[m] - G * bg.M_r[m] / bg.r[m] ** 2) / bg.g[m]) < 1e-12
    # Hydrostatic equilibrium, away from the surface where P falls off a cliff.
    i = bg.x < 0.95
    dPdr = np.gradient(bg.P, bg.r)
    assert np.median(np.abs(dPdr[i] + bg.rho[i] * bg.g[i]) / np.abs(dPdr[i])[i.sum() // 2]) < 1e-3


def test_normalisation(poly3):
    """Test 1 on a second structure, and on per-mode grids."""
    _, efs = poly3
    worst = max(
        abs(2 * ef.omega**2 * ef.inertia_integral() / ef.bg.E_star - 1) for ef in efs.values()
    )
    assert worst < 1e-12


def test_divergence_identity(poly3):
    """Test 2: div.xi vs r^-2 d(r^2 xi_r)/dr - Lambda^2 xi_h/r."""
    _, efs = poly3
    from scipy.interpolate import CubicSpline

    for ef in efs.values():
        b = ef.bg
        m = (b.x > 0.02) & (b.x < 0.95)  # same window as the MESA test
        r = np.where(b.r > 0, b.r, 1.0)
        # Spline derivative, not np.gradient: this grid is 5x coarser than the
        # MESA one near the centre, and a second-order difference is the
        # limiting error there rather than anything in modes.py.
        d = CubicSpline(b.r, b.r**2 * ef.xi_r).derivative()(b.r)
        num = d / r**2 - ef.Lambda2 * ef.xi_h / r
        scale = np.max(np.abs(ef.div_xi[m]))
        assert np.max(np.abs(num[m] - ef.div_xi[m])) / scale < 1e-4


def test_modes_do_not_share_a_grid(poly3):
    """The premise of the union-grid path: GYRE refines per mode here, unlike
    the delta Sct model whose 7137-point grid it never needed to touch."""
    _, efs = poly3
    assert len({id(ef.bg) for ef in efs.values()}) > 1


@pytest.mark.parametrize("keys", CASES3)
def test_permutation_symmetry(poly3, keys):
    """Test 5. On MIXED_L3 the three modes sit on different grids, so this also
    checks that the union grid and the interpolation onto it are permutation
    invariant -- nothing else exercises that."""
    _, efs = poly3
    a, b, c = _efs(efs, keys)
    k = [
        kappa_abc(a, b, c, (0, 0, 0)).kappa,
        kappa_abc(b, c, a, (0, 0, 0)).kappa,
        kappa_abc(c, a, b, (0, 0, 0)).kappa,
    ]
    assert max(abs(x - k[0]) for x in k) / abs(k[0]) < 1e-6


def test_union_grid_is_used_when_grids_differ(poly3):
    _, efs = poly3
    a, b, c = _efs(efs, MIXED_L3)
    assert len({id(a.bg), id(b.bg), id(c.bg)}) > 1
    from coupling.kappa import _on_grid

    grid, _ = _on_grid((a, b, c), 1)
    n_union = len(np.unique(np.concatenate([a.x, b.x, c.x])))
    assert len(grid.r) == n_union
    assert len(grid.r) >= max(len(a.x), len(b.x), len(c.x))


@pytest.mark.parametrize("keys", CASES3)
def test_grid_convergence(poly3, keys):
    _, efs = poly3
    k = grid_convergence(*_efs(efs, keys), (0, 0, 0), refines=(8, 16, 32))
    assert abs(k[32] - k[8]) / abs(k[32]) < 1e-3
    assert abs(k[32] - k[16]) < abs(k[16] - k[8])


@pytest.mark.parametrize("keys", CASES3)
def test_adaptive_result_is_accurate(poly3, keys):
    """The refinement ladder has to reach the converged value on its own."""
    _, efs = poly3
    a, b, c = _efs(efs, keys)
    auto = kappa_abc(a, b, c, (0, 0, 0))
    truth = kappa_abc(a, b, c, (0, 0, 0), refine=32).kappa
    assert auto.converged
    assert abs(auto.kappa - truth) / abs(truth) < 1e-3


def test_g_mode_group_dominance(poly3):
    """Test 10 on a second structure."""
    _, efs = poly3
    r = kappa_abc(*_efs(efs, G3), (0, 0, 0))
    assert max(GROUPS, key=lambda k: abs(r.groups[k])) in ("A56", "A58")
    assert r.outer_fraction(0.8) < 0.1


# ---------------------------------------------------------------- n = 0

KELVIN = {2: 0.8, 3: 2 * 3 * 2 / 7}  # omega^2 = 2l(l-1)/(2l+1)


def test_constant_density(poly0):
    """rho is constant to 2e-15, but dlnrho/dlnr is only zero to ~0.13 at the
    outermost points. It comes from -V_2 x^2/Gamma_1 - As, and at x = 0.998
    those two terms are each ~326: GYRE's As and V_2 disagree at the 4e-4 level
    across the hard n = 0 surface, which is a real inconsistency in the model,
    not float64 roundoff (that would give 3e-14)."""
    bg, _ = poly0
    assert np.ptp(bg.rho) / bg.rho[0] < 1e-9
    d = np.abs(bg.dlnrho_dlnr)
    assert np.median(d) < 1e-14
    assert np.max(d[bg.x < 0.99]) < 1e-6
    assert np.max(d) < 1.0


@pytest.mark.parametrize("l", [2, 3])
def test_kelvin_f_mode_frequency(poly0, l):
    """omega^2 = 2l(l-1)/(2l+1) exactly: the f-mode is incompressible, so
    Gamma_1 never enters and the uniform-density result is not an approximation."""
    _, efs = poly0
    ef = efs[(l, 0)]
    assert ef.omega_dimless**2 == pytest.approx(KELVIN[l], rel=2e-3)


@pytest.mark.parametrize("l", [2, 3])
def test_kelvin_f_mode_eigenfunction(poly0, l):
    """xi_r ~ r^(l-1), xi_h = xi_r/l, div.xi = 0."""
    _, efs = poly0
    ef = efs[(l, 0)]
    b = ef.bg
    m = (b.x > 0.05) & (b.x < 0.95)
    shape = b.x[m] ** (l - 1)
    assert np.ptp(ef.xi_r[m] / shape) / abs(np.mean(ef.xi_r[m] / shape)) < 1e-3
    assert np.max(np.abs(ef.xi_h[m] * l - ef.xi_r[m])) / np.max(np.abs(ef.xi_r[m])) < 1e-3
    assert np.max(np.abs(ef.div_xi[m])) / np.max(np.abs(ef.xi_r[m]) / b.R) < 1e-3


@pytest.mark.parametrize("keys", [[(2, 3), (2, 1), (2, 4)], [(2, 1), (2, 1), (2, 2)]])
def test_A57_is_negligible(poly0, keys):
    """A57 carries the only factor of dlnrho/dlnr, so constant density should
    switch it off. It comes out at 1e-4 of the largest other group rather than
    at zero -- entirely from the surface inconsistency in test_constant_density,
    amplified by xi_r^3. Nothing else tests the A57 code path this way."""
    _, efs = poly0
    r = kappa_abc(*_efs(efs, keys), (0, 0, 0))
    other = max(abs(v) for k, v in r.groups.items() if k != "A57")
    assert abs(r.groups["A57"]) < 1e-3 * other


def test_f_mode_triplet_isolates_the_gravity_groups(poly0):
    """div.xi = 0 kills A55, A56, A58 and, with dlnrho/dlnr = 0, A62 too. An
    (f,f,f) triplet is A59 + A60 + A61 alone."""
    _, efs = poly0
    r = kappa_abc(efs[(2, 0)], efs[(2, 0)], efs[(2, 0)], (0, 0, 0))
    live = {"A59", "A60", "A61"}
    scale = max(abs(r.groups[k]) for k in live)
    assert scale > 0
    for name in GROUPS:
        if name not in live:
            assert abs(r.groups[name]) < 1e-3 * scale
    assert r.kappa == pytest.approx(sum(r.groups[k] for k in live), rel=1e-2)
