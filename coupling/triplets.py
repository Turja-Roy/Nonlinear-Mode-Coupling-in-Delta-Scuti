"""Resonant triplet enumeration.

omega is m-degenerate for a nonrotating star, so the detuning depends only on
(l, n_pg). Enumeration therefore works at that level and attaches the m
combinations afterwards; kappa.py exploits the same split, computing each
radial integral once per RadialTriplet and contracting it with the angular
factors of every m combination.

Sign classes: with all stored omega positive, (+,+,+) can never resonate here
(the lowest three frequencies already sum to 3.3 c/d against a 0.43 c/d cut),
and any two-minus assignment is the negation of a one-minus one, so |Delta| is
unchanged. One mode carries s = -1 — the sum mode, always the highest-frequency
member — and there are three distinct choices per unordered triple.

Slots a, b, c here are bookkeeping, not roles: a is simply the sum mode. Which
modes are parents (gamma < 0, kappa-driven) and which are daughters is decided
by linear stability alone — see `observables.channel`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .angular import l_multisets, m_combos, satisfies_selection_rules
from .modes import Eigenfunction, Mode

DETUNING_CUT_DIMLESS = 0.15  # in units of sqrt(GM/R^3)

Key = tuple[int, int]  # (l, n_pg)


@dataclass(frozen=True)
class RadialTriplet:
    """One (l, n_pg) triple with a fixed sign assignment, m not yet chosen."""

    sum_mode: Key
    pair: tuple[Key, Key]
    omega: tuple[float, float, float]  # signed, rad/s; sum mode first and negative
    delta: float  # rad/s

    @property
    def keys(self) -> tuple[Key, Key, Key]:
        return (self.sum_mode, *self.pair)

    @property
    def ls(self) -> tuple[int, int, int]:
        return tuple(k[0] for k in self.keys)

    @property
    def ns(self) -> tuple[int, int, int]:
        return tuple(k[1] for k in self.keys)

    @property
    def signs(self) -> tuple[int, int, int]:
        return (-1, +1, +1)

    def m_combinations(self) -> list[tuple[int, int, int]]:
        return m_combos(*self.ls)

    def modes(self, efs: dict[Key, Eigenfunction], ms: tuple[int, int, int]) -> tuple[Mode, ...]:
        return tuple(
            Mode(n_pg=k[1], l=k[0], m=m, s=s, ef=efs[k])
            for k, m, s in zip(self.keys, ms, self.signs)
        )

    def __repr__(self) -> str:
        (la, na), (lb, nb), (lc, nc) = self.keys
        return f"RadialTriplet(-({la},{na:+d}) + ({lb},{nb:+d}) + ({lc},{nc:+d}), Delta={self.delta:+.3e} rad/s)"


def enumerate_triplets(
    efs: dict[Key, Eigenfunction],
    cut: float,
    l_max: int = 3,
    sum_keys: set[Key] | None = None,
) -> list[RadialTriplet]:
    """Selection rules first, then |Delta| < cut. `cut` in rad/s.

    `sum_keys` restricts slot a, the sum mode, to those (l, n_pg). The pair
    still ranges over the whole net. That asymmetry is what makes a wide-net
    scan tractable: only a self-excited mode can pump a parametric pair, and the
    driven set is a few dozen modes against a few thousand.
    """
    by_l: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    sum_l: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for l in range(l_max + 1):
        keys = sorted((k for k in efs if k[0] == l), key=lambda k: efs[k].omega)
        if keys:
            by_l[l] = (np.array([k[1] for k in keys]), np.array([efs[k].omega for k in keys]))
            sub = keys if sum_keys is None else [k for k in keys if k in sum_keys]
            if sub:
                sum_l[l] = (np.array([k[1] for k in sub]),
                            np.array([efs[k].omega for k in sub]))

    out: list[RadialTriplet] = []
    for l_multiset in l_multisets(l_max):
        if max(l_multiset) > l_max:
            continue
        for l_p, l_q, l_r in _sum_slot_assignments(l_multiset):
            if l_p not in sum_l or l_q not in by_l or l_r not in by_l:
                continue
            out.extend(_match(efs, sum_l[l_p], by_l, l_p, l_q, l_r, cut))
    return out


def _sum_slot_assignments(l_multiset: tuple[int, int, int]) -> list[tuple[int, int, int]]:
    """Distinct (sum-slot l, l, l); the pair is unordered."""
    seen: set[tuple[int, tuple[int, int]]] = set()
    out = []
    for i in range(3):
        rest = tuple(sorted(l_multiset[:i] + l_multiset[i + 1 :]))
        tag = (l_multiset[i], rest)
        if tag not in seen:
            seen.add(tag)
            out.append((l_multiset[i], *rest))
    return out


def _match(efs, p_modes, by_l, l_p, l_q, l_r, cut) -> list[RadialTriplet]:
    # The sum slot comes in separately: restricting it must not also restrict a
    # pair slot that happens to share its l.
    n_p, w_p = p_modes
    n_q, w_q = by_l[l_q]
    n_r, w_r = by_l[l_r]

    # The pair is unordered, so restrict to the upper triangle when its two
    # members share an l. q == r is kept: that is the self-coupled case.
    iq, ir = np.meshgrid(np.arange(len(n_q)), np.arange(len(n_r)), indexing="ij")
    keep = iq <= ir if l_q == l_r else np.ones_like(iq, dtype=bool)
    iq, ir = iq[keep], ir[keep]
    sums = w_q[iq] + w_r[ir]

    lo = np.searchsorted(w_p, sums - cut, side="left")
    hi = np.searchsorted(w_p, sums + cut, side="right")

    out = []
    for j, (a, b) in enumerate(zip(lo, hi)):
        for ip in range(a, b):
            out.append(
                RadialTriplet(
                    sum_mode=(l_p, int(n_p[ip])),
                    pair=((l_q, int(n_q[iq[j]])), (l_r, int(n_r[ir[j]]))),
                    omega=(-float(w_p[ip]), float(w_q[iq[j]]), float(w_r[ir[j]])),
                    delta=float(sums[j] - w_p[ip]),
                )
            )
    return out


def count_with_m(triplets: list[RadialTriplet]) -> int:
    return sum(len(t.m_combinations()) for t in triplets)


def to_frame(triplets: list[RadialTriplet], with_m: bool = True):
    """Flat table, sum mode in the `a` slot. Persist this before running kappa.

    `with_m` off drops the n_m column. Counting m combinations goes through
    sympy's wigner_3j, which is affordable at l <= 3 and not at l ~ 25.
    """
    import pandas as pd

    return pd.DataFrame(
        [
            {
                "l_a": t.ls[0], "n_a": t.ns[0],
                "l_b": t.ls[1], "n_b": t.ns[1],
                "l_c": t.ls[2], "n_c": t.ns[2],
                "omega_a": t.omega[0], "omega_b": t.omega[1], "omega_c": t.omega[2],
                "delta": t.delta,
                "abs_delta_over_omega_a": abs(t.delta / t.omega[0]),
                **({"n_m": len(t.m_combinations())} if with_m else {}),
            }
            for t in triplets
        ]
    )


def from_frame(df) -> list[RadialTriplet]:
    return [
        RadialTriplet(
            sum_mode=(int(r.l_a), int(r.n_a)),
            pair=((int(r.l_b), int(r.n_b)), (int(r.l_c), int(r.n_c))),
            omega=(float(r.omega_a), float(r.omega_b), float(r.omega_c)),
            delta=float(r.delta),
        )
        for r in df.itertuples()
    ]


def summarise(triplets: list[RadialTriplet], omega_dyn: float) -> dict:
    if not triplets:
        return {"n_radial": 0, "n_with_m": 0}
    delta = np.array([t.delta for t in triplets])
    omega_p = np.array([-t.omega[0] for t in triplets])
    by_ls: dict[tuple[int, int, int], int] = {}
    for t in triplets:
        key = tuple(sorted(t.ls))
        by_ls[key] = by_ls.get(key, 0) + 1
    return {
        "n_radial": len(triplets),
        "n_with_m": count_with_m(triplets),
        "min_abs_delta_over_omega": float(np.min(np.abs(delta) / omega_p)),
        "median_abs_delta_over_omega": float(np.median(np.abs(delta) / omega_p)),
        "l_combinations": dict(sorted(by_ls.items())),
    }


__all__ = [
    "DETUNING_CUT_DIMLESS",
    "RadialTriplet",
    "count_with_m",
    "enumerate_triplets",
    "from_frame",
    "satisfies_selection_rules",
    "summarise",
    "to_frame",
]
