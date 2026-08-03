"""What a triplet does: coupling strength, growth rates, thresholds.

A triplet supports two distinct readings, and they use different modes:

- **direct / combination frequency** — two modes beat and drive the third at
  their sum frequency, q_x = mu q_y q_z. `mu` belongs to the *driven* mode, so
  all three are reported; ordinary combination frequencies sit near mu ~ 4.
- **parametric** — the parent decays into the two daughters once its energy
  exceeds `threshold_energy`.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from .background import Background
from .kappa import KappaResult, kappa_all_m
from .modes import Eigenfunction
from .triplets import RadialTriplet

NONADIABATIC_CYCLES = 2.0 * np.pi


def mu(kappa: float, omega: float, delta: float, gamma: float) -> float:
    """|omega_x kappa| / sqrt(Delta^2 + gamma_x^2) for the driven mode x."""
    return abs(omega * kappa) / np.hypot(delta, gamma)


def parametric_growth_rate(kappa: float, omega_b: float, omega_c: float, q_a: float) -> float:
    """Gamma = 2 kappa sqrt(omega_b omega_c) |q_a|, the daughters' growth rate."""
    return 2.0 * abs(kappa) * np.sqrt(abs(omega_b * omega_c)) * abs(q_a)


def threshold_energy(
    kappa: float,
    omega_b: float,
    omega_c: float,
    gamma_b: float,
    gamma_c: float,
    delta: float,
    E_star: float,
) -> float:
    """Parent energy at which the daughter pair goes unstable.

    Returns 0 if either daughter is driven rather than damped: the pair is
    then linearly unstable and the parametric threshold does not apply.
    """
    if gamma_b <= 0.0 or gamma_c <= 0.0 or kappa == 0.0:
        return 0.0
    detune = 1.0 + delta**2 / (gamma_b + gamma_c) ** 2
    return (gamma_b * gamma_c) / (4.0 * kappa**2 * abs(omega_b * omega_c)) * detune * E_star


def threshold_energy_ceiling(
    kappa: float, omega_b: float, omega_c: float, delta: float, E_star: float
) -> float:
    """Largest `threshold_energy` any gamma pair can give at |Delta| >> gamma.

    In that limit E_th = [gamma_b gamma_c/(gamma_b+gamma_c)^2] Delta^2/(4 k^2 w_b w_c),
    and the bracket is at most 1/4, reached only when the daughters damp at the
    same rate. Unequal damping lowers the threshold, and one undamped daughter
    sends it to zero. So this is an upper bound, which is what to quote while
    gamma is the least certain input.
    """
    if kappa == 0.0:
        return np.inf
    return delta**2 / (16.0 * kappa**2 * abs(omega_b * omega_c)) * E_star


def nonadiabatic_shell(bg: Background, omega: float) -> np.ndarray:
    """Mask where omega t_thermal < 2 pi, i.e. where an adiabatic eigenfunction
    is no longer justified. Sun et al. put this at r/R >~ 0.96 for a delta Sct."""
    return abs(omega) * bg.t_thermal < NONADIABATIC_CYCLES


def nonadiabatic_fraction(res: KappaResult, bg: Background, omega: float) -> float:
    """Share of |dkappa/dr| accumulated inside that shell. A large value is a
    physics caveat on the triplet, not a numerical one."""
    t_th = np.interp(res.r, bg.r, bg.t_thermal)
    w = np.abs(res.dkappa_dr)
    total = np.trapezoid(w, res.r)
    if total == 0.0:
        return 0.0
    m = abs(omega) * t_th < NONADIABATIC_CYCLES
    return float(np.trapezoid(np.where(m, w, 0.0), res.r) / total)


@dataclasses.dataclass(frozen=True)
class TripletObservables:
    ms: tuple[int, int, int]
    kappa: float
    mu: tuple[float, float, float]  # per driven mode
    E_threshold: float  # in erg
    refine: int

    @property
    def mu_max(self) -> float:
        return max(self.mu)


def observables(
    t: RadialTriplet,
    efs: dict[tuple[int, int], Eigenfunction],
    bg: Background,
) -> list[TripletObservables]:
    ks, refine = kappa_all_m(t, efs)
    gam = [efs[k].gamma for k in t.keys]
    out = []
    for ms, k in ks.items():
        out.append(
            TripletObservables(
                ms=ms,
                kappa=k,
                mu=tuple(mu(k, t.omega[i], t.delta, gam[i]) for i in range(3)),
                E_threshold=threshold_energy(
                    k, t.omega[1], t.omega[2], gam[1], gam[2], t.delta, bg.E_star
                ),
                refine=refine,
            )
        )
    return out


def to_frame(
    triplets: list[RadialTriplet],
    efs: dict[tuple[int, int], Eigenfunction],
    bg: Background,
):
    """Ranked table: one row per (triplet, m combination)."""
    import pandas as pd

    rows = []
    for t in triplets:
        gam = [efs[k].gamma for k in t.keys]
        for o in observables(t, efs, bg):
            rows.append(
                {
                    "l_a": t.ls[0], "n_a": t.ns[0], "m_a": o.ms[0],
                    "l_b": t.ls[1], "n_b": t.ns[1], "m_b": o.ms[1],
                    "l_c": t.ls[2], "n_c": t.ns[2], "m_c": o.ms[2],
                    "omega_a": t.omega[0], "omega_b": t.omega[1], "omega_c": t.omega[2],
                    "gamma_a": gam[0], "gamma_b": gam[1], "gamma_c": gam[2],
                    "delta": t.delta,
                    "kappa": o.kappa,
                    "mu_a": o.mu[0], "mu_b": o.mu[1], "mu_c": o.mu[2],
                    "mu_max": o.mu_max,
                    "E_th_over_E_star": o.E_threshold / bg.E_star,
                    "E_th_ceiling_over_E_star": threshold_energy_ceiling(
                        o.kappa, t.omega[1], t.omega[2], t.delta, 1.0
                    ),
                    "detuning_dominated": abs(t.delta) > abs(gam[0]),
                    "refine": o.refine,
                }
            )
    df = pd.DataFrame(rows)
    return df.sort_values("mu_max", ascending=False, ignore_index=True)


__all__ = [
    "NONADIABATIC_CYCLES",
    "TripletObservables",
    "mu",
    "nonadiabatic_fraction",
    "nonadiabatic_shell",
    "observables",
    "parametric_growth_rate",
    "threshold_energy",
    "threshold_energy_ceiling",
    "to_frame",
]
