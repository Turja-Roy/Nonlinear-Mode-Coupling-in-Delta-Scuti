"""What a triplet does: coupling strength, growth rates, thresholds.

Roles follow MW23: a **parent** is self-excited (gamma < 0, kappa-driven), a
**daughter** is damped (gamma > 0) and gets excited nonlinearly. The slots a, b,
c of a RadialTriplet carry no role — a is just the sum mode.

A triplet supports two distinct readings, and they use different modes:

- **direct / combination frequency** — two parents beat and drive one daughter
  at their sum or difference frequency, q_x = mu q_y q_z. `mu` belongs to the
  *daughter*, so all three slots are reported; ordinary combination frequencies
  sit near mu ~ 4.
- **parametric** — one parent decays into two daughters once its energy exceeds
  `threshold_energy`.
"""

from __future__ import annotations

from dataclasses import dataclass

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


CHANNELS = ("parametric", "direct-sum", "direct-diff", "all-driven", "inactive", "all-damped")


def channel(gamma_a: float, gamma_b: float, gamma_c: float) -> str:
    """Which channel the triplet can run, from linear stability alone.

    gamma < 0 is a parent, gamma > 0 a potential daughter. Slot a is the sum
    mode; it is a parent in some channels and the daughter in others. A direct
    channel needs two parents, the parametric one needs a single parent that is
    also the sum mode:

    - parametric: a is the only parent, b and c are its two daughters
    - direct-sum: b and c are parents, they beat and drive daughter a at
      omega_b + omega_c
    - direct-diff: a and one of b, c are parents, driving the remaining
      daughter at the difference frequency
    - all-driven: three parents, no daughter to absorb the energy
    - all-damped: no parent
    - inactive: a lone parent among b, c; its daughters are not in this triplet
    """
    drv = (gamma_a < 0.0, gamma_b < 0.0, gamma_c < 0.0)
    n = sum(drv)
    if n == 3:
        return "all-driven"
    if n == 2:
        return "direct-diff" if drv[0] else "direct-sum"
    if n == 0:
        return "all-damped"
    return "parametric" if drv[0] else "inactive"


def classify_frame(df):
    """`channel` over a ranked table; gamma is m-degenerate, so this is a
    property of the radial triplet."""
    import pandas as pd

    a, b, c = (df[f"gamma_{s}"].to_numpy() < 0.0 for s in "abc")
    n = a.astype(int) + b + c
    out = np.where(n == 3, "all-driven", np.where(n == 0, "all-damped", "inactive"))
    out = np.where((n == 2) & a, "direct-diff", out)
    out = np.where((n == 2) & ~a, "direct-sum", out)
    out = np.where((n == 1) & a, "parametric", out)
    return pd.Series(out, index=df.index, name="channel")


def daughter_index(df) -> np.ndarray:
    """Slot (0, 1, 2) of the daughter a direct channel drives; -1 where no
    direct channel runs, i.e. where the triplet has no single daughter --
    parametric has two, and the rest have none. The mu-argmax slot answers a
    different question -- which mode responds most strongly -- and the two
    disagree on ~13% of rows."""
    b_damped = df["gamma_b"].to_numpy() > 0.0
    cls = classify_frame(df).to_numpy()
    return np.where(
        cls == "direct-sum", 0, np.where(cls == "direct-diff", np.where(b_damped, 1, 2), -1)
    )


@dataclass(frozen=True)
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
    "CHANNELS",
    "NONADIABATIC_CYCLES",
    "TripletObservables",
    "channel",
    "classify_frame",
    "daughter_index",
    "mu",
    "nonadiabatic_fraction",
    "nonadiabatic_shell",
    "observables",
    "parametric_growth_rate",
    "threshold_energy",
    "threshold_energy_ceiling",
    "to_frame",
]
