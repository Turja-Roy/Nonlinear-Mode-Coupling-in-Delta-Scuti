"""Validation suite items 3 and 4."""

import itertools

import numpy as np
import pytest
from scipy.special import roots_legendre, sph_harm_y

from coupling.angular import (
    T,
    angular_factors,
    l_multisets,
    m_combos,
    satisfies_selection_rules,
)


def gaunt_numerical(l_a, l_b, l_c, m_a, m_b, m_c, n=400):
    """int Y_a Y_b Y_c dOmega by quadrature. Sum m = 0 makes the phi integral 2 pi."""
    if m_a + m_b + m_c != 0:
        return 0.0
    cos_t, w = roots_legendre(n)
    theta = np.arccos(cos_t)
    y = np.ones_like(theta, dtype=complex)
    for l, m in ((l_a, m_a), (l_b, m_b), (l_c, m_c)):
        y *= sph_harm_y(l, m, theta, np.zeros_like(theta))
    return float(np.real(2 * np.pi * np.sum(w * y)))


@pytest.mark.parametrize(
    "ls,expected",
    [((0, 0, 0), True), ((1, 1, 1), False), ((0, 0, 2), False), ((1, 1, 2), True), ((2, 2, 2), True)],
)
def test_selection_rules(ls, expected):
    """Test 4."""
    assert satisfies_selection_rules(*ls) is expected


def test_eight_l_combinations():
    combos = l_multisets(3)
    assert len(combos) == 8
    assert combos == [(0, 0, 0), (0, 1, 1), (0, 2, 2), (0, 3, 3), (1, 1, 2), (1, 2, 3), (2, 2, 2), (2, 3, 3)]


def test_gaunt_against_quadrature():
    """Test 3: T from Wigner-3j vs a direct numerical int Y_a Y_b Y_c dOmega."""
    worst = 0.0
    for l_a, l_b, l_c in l_multisets(3):
        for m_a, m_b, m_c in m_combos(l_a, l_b, l_c):
            exact = T(l_a, l_b, l_c, m_a, m_b, m_c)
            num = gaunt_numerical(l_a, l_b, l_c, m_a, m_b, m_c)
            worst = max(worst, abs(exact - num))
    assert worst < 1e-12, f"worst absolute difference {worst:.3e}"


def test_T_vanishes_off_selection_rules():
    assert T(1, 1, 1, 0, 0, 0) == 0.0
    assert T(0, 0, 2, 0, 0, 0) == 0.0
    assert T(1, 1, 2, 1, 0, 0) == 0.0  # sum m != 0


def test_T_permutation_symmetric():
    for l_a, l_b, l_c in l_multisets(3):
        for m_a, m_b, m_c in m_combos(l_a, l_b, l_c):
            ref = T(l_a, l_b, l_c, m_a, m_b, m_c)
            for p in itertools.permutations([(l_a, m_a), (l_b, m_b), (l_c, m_c)]):
                (la, ma), (lb, mb), (lc, mc) = p
                assert T(la, lb, lc, ma, mb, mc) == pytest.approx(ref, abs=1e-14)


def test_F_G_follow_T_permutations():
    """F_a and G_a must track which slot the mode sits in, not its l alone."""
    f = angular_factors(1, 2, 3, 0, 0, 0)
    g = angular_factors(2, 3, 1, 0, 0, 0)
    assert f.F_a == pytest.approx(g.F_c)
    assert f.G_a == pytest.approx(g.G_c)
    assert f.S == pytest.approx(g.S)


def test_radial_triplet_is_trivial():
    """l = 0 everywhere: T = 1/sqrt(4 pi), and every F, G vanishes."""
    f = angular_factors(0, 0, 0, 0, 0, 0)
    assert f.T == pytest.approx(1 / np.sqrt(4 * np.pi))
    assert (f.F_a, f.F_b, f.F_c, f.G_a, f.G_b, f.G_c, f.S) == (0,) * 7
