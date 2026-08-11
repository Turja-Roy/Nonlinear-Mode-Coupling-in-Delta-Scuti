"""Angular integrals for the three-mode coupling coefficient.

T is the Gaunt integral int Y_a Y_b Y_c dOmega for unit-normalised complex
spherical harmonics, int |Y_lm|^2 dOmega = 1 — the same convention the
renormalisation in modes.py assumes. No conjugates: the T3 index convention
puts both frequency signs in the mode list instead.
"""

from __future__ import annotations

from dataclasses import dataclass
import functools
import itertools
import math

from sympy.physics.wigner import wigner_3j


@functools.cache
def wigner3j(l1: int, l2: int, l3: int, m1: int, m2: int, m3: int) -> float:
    return float(wigner_3j(l1, l2, l3, m1, m2, m3))


def satisfies_selection_rules(l_a: int, l_b: int, l_c: int) -> bool:
    """Triangle inequality and even parity; both come from 3j(l_a l_b l_c; 0 0 0)."""
    return (
        abs(l_b - l_c) <= l_a <= l_b + l_c
        and (l_a + l_b + l_c) % 2 == 0
    )


@functools.cache
def T(l_a: int, l_b: int, l_c: int, m_a: int, m_b: int, m_c: int) -> float:
    if m_a + m_b + m_c != 0 or not satisfies_selection_rules(l_a, l_b, l_c):
        return 0.0
    norm = math.sqrt((2 * l_a + 1) * (2 * l_b + 1) * (2 * l_c + 1) / (4 * math.pi))
    return (
        norm
        * wigner3j(l_a, l_b, l_c, m_a, m_b, m_c)
        * wigner3j(l_a, l_b, l_c, 0, 0, 0)
    )


@dataclass(frozen=True)
class AngularFactors:
    T: float
    F_a: float
    F_b: float
    F_c: float
    G_a: float
    G_b: float
    G_c: float
    S: float


@functools.cache
def angular_factors(
    l_a: int, l_b: int, l_c: int, m_a: int, m_b: int, m_c: int
) -> AngularFactors:
    t = T(l_a, l_b, l_c, m_a, m_b, m_c)
    La, Lb, Lc = (l * (l + 1) for l in (l_a, l_b, l_c))

    F_a = 0.5 * t * (Lb + Lc - La)
    F_b = 0.5 * t * (Lc + La - Lb)
    F_c = 0.5 * t * (La + Lb - Lc)
    G_a = 0.25 * t * (La**2 - (Lb - Lc) ** 2)
    G_b = 0.25 * t * (Lb**2 - (Lc - La) ** 2)
    G_c = 0.25 * t * (Lc**2 - (La - Lb) ** 2)

    return AngularFactors(
        T=t,
        F_a=F_a,
        F_b=F_b,
        F_c=F_c,
        G_a=G_a,
        G_b=G_b,
        G_c=G_c,
        S=0.5 * (La * F_a + Lb * F_b + Lc * F_c),
    )


def l_multisets(l_max: int) -> list[tuple[int, int, int]]:
    """Unordered {l_a, l_b, l_c} passing the selection rules. 8 of them for l <= 3."""
    return [
        combo
        for combo in itertools.combinations_with_replacement(range(l_max + 1), 3)
        if satisfies_selection_rules(*combo)
    ]


def m_combos(l_a: int, l_b: int, l_c: int) -> list[tuple[int, int, int]]:
    """Every (m_a, m_b, m_c) with sum zero and T != 0."""
    return [
        (m_a, m_b, -m_a - m_b)
        for m_a in range(-l_a, l_a + 1)
        for m_b in range(-l_b, l_b + 1)
        if abs(m_a + m_b) <= l_c and T(l_a, l_b, l_c, m_a, m_b, -m_a - m_b) != 0.0
    ]
