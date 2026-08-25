"""Amplitude equations for mode networks of any size.

Modes carry complex amplitudes q_x and unsigned frequencies. The normalisation
2 omega^2 int rho |xi|^2 d3x = E_star is chosen precisely so that

    E_x = |q_x|^2 E_star,

with no frequency factor; N_x = E_x / omega_x is the Manley-Rowe action.

One Coupling holds one triplet, with `sum_slot` naming the mode at
omega_sum ~= omega_1 + omega_2. The role a mode plays in a triplet, not the
mode itself, fixes the conjugations:

    sum slot:  qdot_x + (i w_x + g_x) q_x = i s_x w_x kappa  q_u q_v
    pair slot: qdot_x + (i w_x + g_x) q_x = i s_x w_x kappa* q_sum q_other*

with the combinatorial factor s_x counting the orderings of the two partners:
2 when they are distinct modes, 1 when they are the same mode. That is what
makes the triplet conserve energy on resonance. MW25's Appendix absorbs the 2
into kappa, so check which convention a quoted kappa is in before comparing.

Carrying the role per coupling rather than as a sign on omega is what lets a
mode be a daughter in one triplet and the pump in another -- the mixed
direct/parametric networks of MW25 sec 5, which no single sign per mode can
express.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class NetworkMode:
    name: str
    omega: float  # rad/s, > 0
    gamma: float  # rad/s, > 0 damping


@dataclass(frozen=True)
class Coupling:
    """One triplet, as indices into the network's mode list."""

    idx: tuple[int, int, int]
    kappa: float
    sum_slot: int = 0

    def terms(self) -> list[tuple[int, float, tuple[int, int], bool]]:
        """(mode, factor, partners, is_sum), once per distinct mode in the triplet.

        For a pair mode the partners come back sum-mode first, which is the
        order the amplitude equation conjugates in. A mode occupying two slots
        still gets a single equation term; the doubling it would pick up from a
        naive slot loop is exactly the doubling already carried by its two
        distinct partners. A repeated mode is always in the pair slots -- a
        mode coincident with its own sum slot would need a zero-frequency
        partner.
        """
        out = []
        for i in dict.fromkeys(self.idx):
            is_sum = i == self.idx[self.sum_slot]
            slot = self.sum_slot if is_sum else self.idx.index(i)
            rest = [j for j in range(3) if j != slot]
            if not is_sum:
                rest.sort(key=lambda j: j != self.sum_slot)
            others = tuple(self.idx[j] for j in rest)
            out.append((i, 1.0 if others[0] == others[1] else 2.0, others, is_sum))
        return out


class Network:
    def __init__(self, modes: list[NetworkMode], couplings: list[Coupling]):
        self.modes = modes
        self.couplings = couplings
        self.names = [m.name for m in modes]
        self.omega = np.array([m.omega for m in modes])
        self.gamma = np.array([m.gamma for m in modes])

    def detuning(self, cp: Coupling) -> float:
        """Delta = (pair frequencies) - (sum-slot frequency)."""
        w = self.omega[list(cp.idx)]
        return float(w.sum() - 2 * w[cp.sum_slot])

    def rhs(self, t: float, q: np.ndarray) -> np.ndarray:
        """Lab frame. Carries the full oscillation, so only useful for short runs."""
        out = -(1j * self.omega + self.gamma) * q
        for cp in self.couplings:
            for i, s, (u, v), is_sum in cp.terms():
                pre = 1j * s * self.omega[i]
                if is_sum:
                    out[i] += pre * cp.kappa * q[u] * q[v]
                else:
                    out[i] += pre * np.conj(cp.kappa) * q[u] * np.conj(q[v])
        return out

    def rhs_rotating(self, t: float, A: np.ndarray) -> np.ndarray:
        """q_x = A_x exp(-i omega_x t). The fast oscillation drops out and the
        coupling keeps only exp(-+ i Delta t) -- autonomous at exact resonance.

        Pulsation periods are hours and amplitude modulation is >= 100 d, so
        this is the only form in which the long integrations are affordable.
        """
        out = -self.gamma * A
        for cp in self.couplings:
            d = self.detuning(cp)
            e_sum, e_pair = np.exp(-1j * d * t), np.exp(1j * d * t)
            for i, s, (u, v), is_sum in cp.terms():
                pre = 1j * s * self.omega[i]
                if is_sum:
                    out[i] += pre * cp.kappa * A[u] * A[v] * e_sum
                else:
                    out[i] += pre * np.conj(cp.kappa) * A[u] * np.conj(A[v]) * e_pair
        return out

    def to_lab(self, t: np.ndarray, A: np.ndarray) -> np.ndarray:
        return A * np.exp(-1j * np.outer(self.omega, t))

    def energy(self, q: np.ndarray) -> np.ndarray:
        """E_x / E_star."""
        return np.abs(q) ** 2

    def action(self, q: np.ndarray) -> np.ndarray:
        """N_x = E_x / omega_x, in units of E_star; positive.

        One triplet moves f quanta out of its sum mode and f into each pair
        mode, so N_sum + N_pair and the difference of the two pair actions are
        the invariants. In a network the sums run over every coupling a mode
        takes part in.
        """
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
        frozen: tuple[int, ...] = (),
    ):
        """Returns (t, A) in the rotating frame, or (t, q) with rotating=False.
        Both give the same |amplitude|, so energies and actions are unaffected.

        `e_max` stops the run when the total energy crosses it, so a runaway
        shows up as t[-1] < t_end. A network whose modes are all linearly
        driven has no bounded state, and without a cap the integrator simply
        grinds against the blow-up.

        `frozen` pins those modes at their initial amplitude: the physical
        picture for linearly driven modes whose saturation lies outside the
        network. Only meaningful in the rotating frame.
        """
        q0 = np.asarray(q0, dtype=complex)
        n = len(q0)
        # Amplitudes are ~1e-6, so a fixed atol would dominate the error budget.
        atol = 1e-14 * float(np.max(np.abs(q0))) if atol is None else atol
        rhs = self.rhs_rotating if rotating else self.rhs
        hold = np.array(frozen, dtype=int)

        def f(t, y):
            z = rhs(t, y[:n] + 1j * y[n:])
            z[hold] = 0.0
            return _split(z)

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
    """Sum mode first, as a RadialTriplet orders it."""
    modes = [NetworkMode(nm, w, g) for nm, w, g in zip(names, omega, gamma)]
    return Network(modes, [Coupling(idx=(0, 1, 2), kappa=kappa, sum_slot=0)])


def self_coupled(
    omega_a: float,
    omega_d: float,
    gamma_a: float,
    gamma_d: float,
    kappa: float,
) -> Network:
    """Sum-slot mode a decaying into two copies of the same mode d."""
    modes = [NetworkMode("a", omega_a, gamma_a), NetworkMode("d", omega_d, gamma_d)]
    return Network(modes, [Coupling(idx=(0, 1, 1), kappa=kappa, sum_slot=0)])


def from_triplets(triplets, efs, ms_list=None, gamma_override=None) -> Network:
    """Network over the union of several RadialTriplets.

    A mode shared between triplets is one node whatever role it plays in each,
    so mixed direct/parametric topologies come through unchanged.
    """
    from .kappa import kappa_all_m

    if ms_list is None:
        ms_list = [None] * len(triplets)
    index: dict[tuple[int, int], int] = {}
    modes: list[NetworkMode] = []
    couplings: list[Coupling] = []

    for t, ms in zip(triplets, ms_list):
        ks, _ = kappa_all_m(t, efs)
        ms = ms if ms is not None else max(ks, key=lambda m: abs(ks[m]))
        idx = []
        for key in t.keys:  # sum mode first
            if key not in index:
                gamma = (gamma_override or {}).get(key, efs[key].gamma)
                index[key] = len(modes)
                modes.append(NetworkMode(f"({key[0]},{key[1]:+d})", efs[key].omega, gamma))
            idx.append(index[key])
        couplings.append(Coupling(idx=tuple(idx), kappa=ks[ms], sum_slot=0))

    return Network(modes, couplings)


__all__ = [
    "Coupling",
    "Network",
    "NetworkMode",
    "from_triplets",
    "self_coupled",
    "three_mode",
]
