#!/usr/bin/env python3
"""Integrate a mixed direct/parametric network to a bounded state.

Direct and parametric coupling act on the same modes (MW25 sec 5). Which mode
carries the parametric leg decides whether the network has a bounded state at
all:

    daughter   a + b -> c,  c -> d + d
    parent     a + b -> c,  a -> d + d,  b -> d + d
    pair       a + b -> c,  a -> d1 + d2,  b -> d1 + d2

Parents are self-excited (gamma < 0, MW23 sec 2) and nothing else in the
`daughter` topology can absorb that flux, so its parents run away and are only
integrable as a fixed-amplitude pump -- the physical picture for modes whose
saturation lies outside the network. In `parent` and `pair` the parametric leg
drains the parents themselves, which is MW25's stabilisation mechanism and the
relaxation-oscillation limit cycle of their Figs. 3-5; those run unpumped.

Topology is data here, not equations: each entry of TOPOLOGIES returns a mode
list and a list of Couplings, and `coupling.network.Network` integrates any
number of them. A five- or ten-mode network is another entry, not new code.

Selection is not on growth rate alone. Gamma/gamma_d ranks the strength of the
instability, but the cost of the integration is set by Delta: the rotating-frame
solution oscillates at Delta while the amplitudes evolve on 1/gamma, and a
system with Gamma << |Delta| needs ~|Delta|/Gamma oscillations per e-folding.
Only near-resonant systems are integrable, so the run picks the best candidate
with Gamma_par comparable to or above |Delta_param|.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from coupling.network import Coupling, Network, NetworkMode
from coupling.observables import mu

CD = 7.27220522e-5
YR = 3.15576e7


def _modes(r, labels):
    """NetworkModes for the (l, n) columns named by `labels`, frequencies in c/d."""
    return [
        NetworkMode(f"{s}({int(r[f'l_{s}'])},{int(r[f'n_{s}']):+d})",
                    float(r[f"f_{s}"]) * CD, float(r[f"gamma_{s}"]))
        for s in labels
    ]


def _direct(r) -> Coupling:
    """a + b -> c. The sum branch puts c in the sum slot, the difference branch
    the higher-frequency parent, which is the assignment four_mode_search made
    when it chose the kappa integral.
    """
    if r.comb == "sum":
        return Coupling(idx=(2, 0, 1), kappa=float(r.kappa_direct), sum_slot=0)
    hi = 0 if r.f_a > r.f_b else 1
    return Coupling(idx=(hi, 1 - hi, 2), kappa=float(r.kappa_direct), sum_slot=0)


def _daughter(r):
    modes = _modes(r, "abcd")
    return modes, [_direct(r), Coupling(idx=(2, 3, 3), kappa=float(r.kappa_param))], (0, 1)


def _parent(r):
    modes = _modes(r, "abcd")
    return modes, [_direct(r),
                   Coupling(idx=(0, 3, 3), kappa=float(r.kappa_param_a)),
                   Coupling(idx=(1, 3, 3), kappa=float(r.kappa_param_b))], ()


def _pair(r):
    modes = _modes(r, "abcde")
    return modes, [_direct(r),
                   Coupling(idx=(0, 3, 4), kappa=float(r.kappa_param_a)),
                   Coupling(idx=(1, 3, 4), kappa=float(r.kappa_param_b))], ()


# name -> (builder, modes frozen by default under --pump)
TOPOLOGIES = {"daughter": _daughter, "parent": _parent, "pair": _pair}


def drive_amplitude(r, topology: str) -> float:
    """Parent amplitude to start from, in units where E = |q|^2 E_star.

    `daughter`: E_c scales as q^4, so the amplitude that puts the direct
    daughter at its parametric threshold follows the quarter power. Driving
    above it automatically puts Gamma above |Delta_p|/2, which is what makes
    the integration affordable at all.

    `parent`/`pair`: the parent crosses its own threshold under kappa-mechanism
    growth, so start below it and let the linear driving carry it across.
    """
    if topology == "daughter":
        return 1e-6 * float(r.E_c_over_E_th) ** -0.25
    return float(r.E_th_par) ** 0.5


def _check_slaving(net: Network, q_thr: float) -> None:
    """Slaving is a sub-threshold statement: above E_th the parametric leg
    drains c and |q_c| sits below the slaved value by construction. Verify it
    where it should hold, at half the threshold amplitude, where d decays and c
    relaxes cleanly. The response is s_c mu q_a q_b with s_c = 2 for distinct
    parents -- the bare mu misses the combinatorial factor of the amplitude
    equation, which is what this measurement caught.
    """
    direct = net.couplings[0]
    d1 = net.detuning(direct)
    slaved = 2.0 * mu(direct.kappa, net.omega[2], d1, net.gamma[2]) * (0.5 * q_thr) ** 2
    q0 = np.full(len(net.modes), 1e-4 * q_thr, dtype=complex)
    q0[[i for i, m in enumerate(net.modes) if m.gamma < 0]] = 0.5 * q_thr
    t_end = max(8.0 / np.hypot(d1, net.gamma[2]), 6.0 / net.gamma[2])
    _, A = net.integrate(q0, t_end, n_out=4000, rtol=1e-10, frozen=(0, 1))
    half = A.shape[1] // 2
    err = np.abs(np.abs(A[2, half:]) - slaved) / slaved
    print(f"\n  q_c vs 2 mu q_a q_b at 0.5x threshold (sub-threshold, slaved): "
          f"median {np.median(err):.3e}")
    print(f"  E_d there: start {abs(A[3, 0])**2:.3e} -> end {abs(A[3, -1])**2:.3e} "
          f"(decays, as it must below threshold)")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--csv", default="out/four_mode_dsct_M2.0.csv")
    p.add_argument("--out", default="out/mixed_network.csv")
    p.add_argument("--topology", default="auto", choices=["auto", *TOPOLOGIES])
    p.add_argument("--q-parent", type=float, default=None,
                   help="parent amplitude; default is set by the topology")
    p.add_argument("--drive", type=float, default=None,
                   help="multiple of the default parent amplitude "
                        "(3 for `daughter`, 0.5 for the self-saturating topologies)")
    p.add_argument("--pump", action="store_true", default=None,
                   help="freeze the parents at their initial amplitude; on by "
                        "default for `daughter`, whose parents cannot saturate")
    p.add_argument("--no-pump", dest="pump", action="store_false")
    p.add_argument("--e-folds", type=float, default=None,
                   help="integration length in e-folds of the slowest rate")
    a = p.parse_args()

    df = pd.read_csv(a.csv)
    topology = a.topology
    if topology == "auto":
        topology = df.branch.iloc[0] if "branch" in df else "daughter"
    if "branch" in df:
        df = df[df.branch == topology]
    if not len(df):
        print(f"no {topology} candidates in {a.csv}")
        return 1

    # Rank on the detuned threshold, not on Gamma/gamma_d: the latter drops the
    # [1 + Delta^2/(2 gamma_d)^2] factor, which is ~(89)^2 at the median.
    rank = "E_c_over_E_th" if topology == "daughter" else "E_par_over_E_th"
    df = df.sort_values(rank, ascending=False)
    print(f"{len(df)} {topology} candidates, {int((df[rank] > 1).sum())} above "
          f"threshold at q = {1e-6:g}")
    r = df.iloc[0]

    modes, couplings, pumped = TOPOLOGIES[topology](r)
    net = Network(modes, couplings)
    print("\nsystem: " + "  ".join(m.name for m in modes))
    for m in modes:
        print(f"  {m.name:<14s} f {m.omega / CD:8.4f} c/d   gamma {m.gamma:+.3e} s^-1")
    for cp in couplings:
        names = "  ".join(modes[i].name for i in cp.idx)
        print(f"  leg [{names}]  kappa {cp.kappa:+.3f}  "
              f"Delta {net.detuning(cp) / CD:+.3e} c/d")
    print(f"  {rank} {r[rank]:.3e} at q = 1e-6")

    pump = (topology == "daughter") if a.pump is None else a.pump
    frozen = pumped if pump else ()
    drive = a.drive if a.drive is not None else (3.0 if topology == "daughter" else 0.5)
    q0_par = a.q_parent if a.q_parent else drive * drive_amplitude(r, topology)
    print(f"\nparent amplitude {q0_par:.3e} ({drive:g}x the default for {topology})")

    # Seeds must sit far enough above atol that a sub-threshold decay does not
    # underflow them: a daughter that reaches zero cannot be re-excited when
    # the parent later crosses, and the parents then overshoot by orders of
    # magnitude before roundoff noise regrows it.
    q0 = np.full(len(modes), 1e-3 * q0_par, dtype=complex)
    q0[[i for i, m in enumerate(modes) if m.gamma < 0]] = q0_par

    rates = [abs(m.gamma) for m in modes]
    if topology == "daughter":
        # The slowest scale is c's off-resonant relaxation, 1/|Delta + i gamma_c|.
        rates = [np.hypot(net.detuning(couplings[0]), modes[2].gamma), modes[3].gamma]
    n_folds = a.e_folds if a.e_folds else (8.0 if topology == "daughter" else 40.0)
    t_end = n_folds / min(rates)
    print(f"{'pumped ' if pump else ''}system, t_end {t_end / YR:.2f} yr"
          + (f", frozen: {', '.join(modes[i].name for i in frozen)}" if frozen else ""))

    # A network with no way to absorb the kappa-mechanism flux has no bounded
    # state; the cap turns that into a finished run instead of a step size
    # shrinking without limit.
    t, A = net.integrate(q0, t_end, n_out=4000, rtol=1e-10, frozen=frozen,
                         e_max=(30.0 * q0_par) ** 2)
    E = net.energy(A)
    if t[-1] < t_end:
        print(f"\n  RUNAWAY: amplitudes passed 30x the drive at "
              f"t = {t[-1] / YR:.3f} yr -- integration stopped")

    print(f"\nintegrated {t[-1] / YR:.3e} yr")
    for m, e in zip(modes, E):
        print(f"  E_{m.name:<14s} start {e[0]:.3e}  end {e[-1]:.3e}  max {e.max():.3e}")

    half = len(t) // 2
    for i, m in enumerate(modes):
        if m.gamma > 0 and E[i][half:].mean() > 0:
            swing = np.ptp(E[i, half:]) / E[i, half:].mean()
            print(f"  {m.name} energy swing over the second half: {swing:.3f}")

    if topology == "daughter" and pump:
        _check_slaving(net, q0_par / drive)

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"t_yr": t / YR, **{f"E_{m.name}": e for m, e in zip(modes, E)}}).to_csv(
        out, index=False)
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
