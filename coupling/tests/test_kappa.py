import numpy as np
import pytest

from coupling.kappa import (
    A56_parts,
    GROUPS,
    grid_convergence,
    kappa_abc,
    kappa_for_triplet,
    prepare,
)

# (l, n_pg) triples spanning the regimes kappa has to survive.
P_TRIPLET = [(2, 5), (2, 15), (2, 6)]
G_TRIPLET = [(1, -13), (1, -11), (2, -10)]
HIGH_N_G = [(1, -18), (1, -19), (2, -20)]
RADIAL = [(0, 9), (0, 1), (0, 6)]
CASES = [P_TRIPLET, G_TRIPLET, HIGH_N_G, RADIAL]


def _efs(efs, keys):
    return [efs[tuple(k)] for k in keys]


def _ms(keys):
    return (0, 0, 0)


@pytest.mark.parametrize("keys", CASES)
def test_permutation_symmetry(efs, keys):
    """Test 5: kappa_abc = kappa_bca = kappa_cab, each evaluated from scratch."""
    a, b, c = _efs(efs, keys)
    m_a, m_b, m_c = _ms(keys)
    k = [
        kappa_abc(a, b, c, (m_a, m_b, m_c)).kappa,
        kappa_abc(b, c, a, (m_b, m_c, m_a)).kappa,
        kappa_abc(c, a, b, (m_c, m_a, m_b)).kappa,
    ]
    assert max(abs(x - k[0]) for x in k) / abs(k[0]) < 1e-3


@pytest.mark.parametrize("keys", CASES)
def test_grid_convergence(efs, keys):
    """Test 6: kappa stable under 2x and 4x refinement."""
    a, b, c = _efs(efs, keys)
    k = grid_convergence(a, b, c, _ms(keys))
    assert abs(k[4] - k[1]) / abs(k[4]) < 1e-3
    assert abs(k[4] - k[2]) < abs(k[2] - k[1])


@pytest.mark.parametrize("keys", CASES)
def test_simpson_agrees(efs, keys):
    a, b, c = _efs(efs, keys)
    assert kappa_abc(a, b, c, _ms(keys)).quadrature_residual < 1e-3


def test_cumulative_shape(efs):
    """Test 7: p-mode triplets build up outside, g-mode triplets near the core."""
    p = kappa_abc(*_efs(efs, P_TRIPLET), _ms(P_TRIPLET))
    g = kappa_abc(*_efs(efs, HIGH_N_G), _ms(HIGH_N_G))
    assert p.outer_fraction(0.8) > 0.8
    assert g.outer_fraction(0.8) < 0.2


@pytest.mark.parametrize("keys", CASES)
def test_no_catastrophic_cancellation(efs, keys):
    """Test 8: max |kappa(<r)| / |kappa(R)| stays small."""
    assert kappa_abc(*_efs(efs, keys), _ms(keys)).cancellation < 1e2


@pytest.mark.parametrize("keys", CASES)
def test_magnitude(efs, keys):
    """Test 9: order unity to ~100, as in MW23."""
    assert 1e-3 < abs(kappa_abc(*_efs(efs, keys), _ms(keys)).kappa) < 1e3


def test_group_dominance(efs):
    """Test 10: A56 leads for high-order g-modes, and inside it the Lambda^2 term."""
    a, b, c = _efs(efs, HIGH_N_G)
    ms = _ms(HIGH_N_G)
    r = kappa_abc(a, b, c, ms)
    assert max(GROUPS, key=lambda k: abs(r.groups[k])) == "A56"

    # Lambda xi_h >> xi_r is a pointwise statement. Both pieces oscillate about
    # zero, so their integrals compare their cancellation, not their size.
    grid, fields, ang = prepare(a, b, c, ms)
    hor, rad = A56_parts(grid, fields, ang)
    cavity = (grid.r > 0.05 * grid.r[-1]) & (grid.r < 0.6 * grid.r[-1])
    assert np.median(np.abs(hor[cavity] / rad[cavity])) > 3.0
    assert np.mean(np.abs(hor[cavity]) > np.abs(rad[cavity])) > 0.75


def test_radial_modes_have_no_horizontal_groups(efs):
    """l = 0 kills xi_h and every F, G, so A59-A61 must vanish identically."""
    r = kappa_abc(*_efs(efs, RADIAL), _ms(RADIAL))
    for name in ("A59", "A60", "A61"):
        assert r.groups[name] == 0.0


def test_cumulative_ends_at_kappa(efs):
    r = kappa_abc(*_efs(efs, G_TRIPLET), _ms(G_TRIPLET))
    assert r.cumulative[-1] == pytest.approx(r.kappa, rel=1e-12)
    assert r.r[-1] == pytest.approx(efs[(1, -13)].bg.R, rel=1e-12)


def test_basis_matches_the_verbatim_groups(efs, bg):
    """The m-factorised path and the literal A55-A62 sum must agree exactly."""
    from coupling.kappa import kappa_all_m
    from coupling.triplets import DETUNING_CUT_DIMLESS, enumerate_triplets

    ts = enumerate_triplets(efs, DETUNING_CUT_DIMLESS * bg.omega_dyn)
    for t in ts[::700]:
        for ms, k in kappa_all_m(t, efs, refine=1)[0].items():
            ref = kappa_for_triplet(t, efs, ms, refine=1).kappa
            assert k == pytest.approx(ref, rel=1e-10, abs=1e-30)


NEAR_NULL = [(1, -5), (1, -17), (2, -15)]  # |kappa| ~ 2.5e-4, L1/|kappa| ~ 4e4


def test_near_null_triplet_is_refined_until_converged(efs):
    """The 0.8% of triplets whose integral, not integrand, cancels. Trapezoid is
    2% off here; the spline rule plus escalation brings it under tol."""
    a, b, c = _efs(efs, NEAR_NULL)
    limit = kappa_abc(a, b, c, (0, 0, 0), refine=32).kappa
    auto = kappa_abc(a, b, c, (0, 0, 0))
    assert auto.refine > 1
    assert abs(auto.kappa - limit) / abs(limit) < 1e-4
    assert abs(np.trapezoid(auto.dkappa_dr, auto.r) - limit) / abs(limit) > 1e-3


def test_sign_assignment_does_not_matter(efs, bg):
    """Only omega^2 appears, so kappa is shared by all three sign classes."""
    from coupling.triplets import DETUNING_CUT_DIMLESS, enumerate_triplets

    t = enumerate_triplets(efs, DETUNING_CUT_DIMLESS * bg.omega_dyn)[0]
    assert all(np.isfinite(kappa_abc(*_efs(efs, t.keys), ms).kappa) for ms in t.m_combinations())
