"""Amplitude equations for small mode networks.

Modes carry complex amplitudes q_x. The normalisation
2 omega^2 int rho |xi|^2 d3x = E_star is chosen precisely so that

    E_x = |q_x|^2 E_star,

with no frequency factor; N_x = E_x / omega_x is the Manley-Rowe action. One
Coupling holds one triplet:

    qdot_x + (i omega_x + gamma_x) q_x = i s_x omega_x kappa* prod_{y != x} q_y*

with omega signed and the combinatorial factor s_x counting the orderings of
the two partners: 2 when they are distinct modes, 1 when they are the same
mode. That is what makes sum_x s_x omega_x = 0 on resonance, which is energy
conservation. MW25's Appendix absorbs the 2 into kappa, so check which
convention a quoted kappa is in before comparing.
"""

from __future__ import annotations

import dataclasses

import numpy as np
from scipy.integrate import solve_ivp


@dataclasses.dataclass(frozen=True)
class NetworkMode:
    name: str
    omega: float  # rad/s, signed
    gamma: float  # rad/s, > 0 damping


@dataclasses.dataclass(frozen=True)
class Coupling:
    """One triplet, as indices into the network's mode list."""

    idx: tuple[int, int, int]
    kappa: float

    def terms(self) -> list[tuple[int, float, tuple[int, int]]]:
        """(mode, factor, partners), once per distinct mode in the triplet.

        A mode occupying two slots still gets a single equation term; the
        doubling it would pick up from a naive slot loop is exactly the
        doubling already carried by its two distinct partners.
        """
        out = []
        for i in dict.fromkeys(self.idx):
            slot = self.idx.index(i)
            others = tuple(self.idx[j] for j in range(3) if j != slot)
            out.append((i, 1.0 if others[0] == others[1] else 2.0, others))
        return out


class Network:
    def __init__(self, modes: list[NetworkMode], couplings: list[Coupling]):
        self.modes = modes
        self.couplings = couplings
        self.omega = np.array([m.omega for m in modes])
        self.gamma = np.array([m.gamma for m in modes])

    def detuning(self, cp: Coupling) -> float:
        return float(self.omega[list(cp.idx)].sum())

    def rhs(self, t: float, q: np.ndarray) -> np.ndarray:
        """Lab frame. Carries the full oscillation, so only useful for short runs."""
        out = -(1j * self.omega + self.gamma) * q
        for cp in self.couplings:
            for i, s, (u, v) in cp.terms():
                out[i] += 1j * s * self.omega[i] * np.conj(cp.kappa) * np.conj(q[u] * q[v])
        return out

    def rhs_rotating(self, t: float, A: np.ndarray) -> np.ndarray:
        """q_x = A_x exp(-i omega_x t). The fast oscillation drops out and the
        coupling keeps only exp(i Delta t) -- autonomous at exact resonance.

        Pulsation periods are hours and amplitude modulation is >= 100 d, so
        this is the only form in which the long integrations are affordable.
        """
        out = -self.gamma * A
        for cp in self.couplings:
            phase = np.exp(1j * self.detuning(cp) * t)
            for i, s, (u, v) in cp.terms():
                out[i] += (
                    1j * s * self.omega[i] * np.conj(cp.kappa) * np.conj(A[u] * A[v]) * phase
                )
        return out

    def to_lab(self, t: np.ndarray, A: np.ndarray) -> np.ndarray:
        return A * np.exp(-1j * np.outer(self.omega, t))

    def energy(self, q: np.ndarray) -> np.ndarray:
        """E_x / E_star."""
        return np.abs(q) ** 2

    def action(self, q: np.ndarray) -> np.ndarray:
        """N_x = E_x / omega_x, in units of E_star; signed with omega."""
        q = np.asarray(q)
        w = self.omega.reshape((-1,) + (1,) * (q.ndim - 1))
        return self.energy(q) / w

    def integrate(
        self,
        q0,
        t_end: float,
        n_out: int = 2000,
        rtol=1e-12,
        atol=None,
        rotating: bool = True,
        e_max: float | None = None,
    ):
        """Returns (t, A) in the rotating frame, or (t, q) with rotating=False.
        Both give the same |amplitude|, so energies and actions are unaffected.

        `e_max` stops the run when the total energy crosses it. A network whose
        modes are all linearly driven has no bounded state, and without a cap
        the integrator simply grinds against the blow-up.
        """
        q0 = np.asarray(q0, dtype=complex)
        n = len(q0)
        # Amplitudes are ~1e-6, so a fixed atol would dominate the error budget.
        atol = 1e-14 * float(np.max(np.abs(q0))) if atol is None else atol
        rhs = self.rhs_rotating if rotating else self.rhs

        def f(t, y):
            return _split(rhs(t, y[:n] + 1j * y[n:]))

        events = None
        if e_max is not None:
            def cap(t, y):
                return float(np.sum(y**2) - e_max)

            cap.terminal = True
            cap.direction = 1.0
            events = [cap]

        sol = solve_ivp(
            f,
            (0.0, t_end),
            _split(q0),
            method="DOP853",
            t_eval=np.linspace(0.0, t_end, n_out),
            rtol=rtol,
            atol=atol,
            events=events,
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.t, sol.y[:n] + 1j * sol.y[n:]


def _split(z: np.ndarray) -> np.ndarray:
    return np.concatenate([z.real, z.imag])


def three_mode(
    omega: tuple[float, float, float],
    gamma: tuple[float, float, float],
    kappa: float,
    names: tuple[str, str, str] = ("a", "b", "c"),
) -> Network:
    modes = [NetworkMode(nm, w, g) for nm, w, g in zip(names, omega, gamma)]
    return Network(modes, [Coupling(idx=(0, 1, 2), kappa=kappa)])


def self_coupled(
    omega_a: float,
    omega_d: float,
    gamma_a: float,
    gamma_d: float,
    kappa: float,
) -> Network:
    """Parent a decaying into two copies of the same daughter d."""
    modes = [NetworkMode("a", omega_a, gamma_a), NetworkMode("d", omega_d, gamma_d)]
    return Network(modes, [Coupling(idx=(0, 1, 1), kappa=kappa)])


def from_triplets(triplets, efs, ms_list=None, gamma_override=None) -> Network:
    """Network over the union of several RadialTriplets.

    Each triplet's parent enters with a negative frequency, so every coupling's
    detuning is the one enumerate_triplets already filtered on. A mode shared
    between triplets keeps whichever sign its first appearance fixed; a triplet
    that would need the opposite sign for a shared mode is rejected.
    """
    from .kappa import kappa_all_m

    if ms_list is None:
        ms_list = [None] * len(triplets)
    index: dict[tuple[int, int], int] = {}
    modes: list[NetworkMode] = []
    signs: dict[tuple[int, int], int] = {}
    couplings: list[Coupling] = []

    for t, ms in zip(triplets, ms_list):
        ks, _ = kappa_all_m(t, efs)
        ms = ms if ms is not None else max(ks, key=lambda m: abs(ks[m]))
        idx = []
        for key, s in zip(t.keys, t.signs):
            if key not in index:
                gamma = (gamma_override or {}).get(key, efs[key].gamma)
                index[key] = len(modes)
                signs[key] = s
                modes.append(NetworkMode(f"({key[0]},{key[1]:+d})", s * efs[key].omega, gamma))
            elif signs[key] != s:
                raise ValueError(f"{key} needs both frequency signs; split the network")
            idx.append(index[key])
        couplings.append(Coupling(idx=tuple(idx), kappa=ks[ms]))

    return Network(modes, couplings)


__all__ = [
    "Coupling",
    "Network",
    "NetworkMode",
    "from_triplets",
    "self_coupled",
    "three_mode",
]
