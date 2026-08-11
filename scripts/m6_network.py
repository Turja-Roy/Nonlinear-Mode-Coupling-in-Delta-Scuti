#!/usr/bin/env python3
"""M6: integrate a four-mode mixed system to a bounded state.

Two legs share the direct daughter c:

    P -> Q + c      (P the higher-frequency driven parent)
    c -> d + d      (self-coupled parametric daughter)

c therefore enters one coupling as a daughter and the other as a parent, which
network.py's one-sign-per-mode Network cannot hold -- it raises on exactly this.
The equations are written out here instead, in the same convention: parent term
carries exp(-i Delta t), daughter terms exp(+i Delta t), and the combinatorial
factor is 2 for distinct partners, 1 for coincident ones.

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
from scipy.integrate import solve_ivp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

CD = 7.27220522e-5
YR = 3.15576e7


def rhs(t, y, w, gam, k1, k2, d1, d2):
    """y = [A_P, A_Q, A_c, A_d] as complex amplitudes."""
    P, Q, c, d = y
    wP, wQ, wc, wd = w
    e1p, e1d = np.exp(-1j * d1 * t), np.exp(1j * d1 * t)
    e2p, e2d = np.exp(-1j * d2 * t), np.exp(1j * d2 * t)
    return np.array([
        -gam[0] * P + 2j * wP * k1 * Q * c * e1p,
        -gam[1] * Q + 2j * wQ * k1 * P * np.conj(c) * e1d,
        -gam[2] * c + 2j * wc * k1 * P * np.conj(Q) * e1d + 1j * wc * k2 * d * d * e2p,
        -gam[3] * d + 2j * wd * k2 * c * np.conj(d) * e2d,
    ])


def integrate(y0, t_end, w, gam, k1, k2, d1, d2, n_out=4000, rtol=1e-10,
              pump=False, amp_cap=None):
    """pump=True freezes P and Q: the parents are treated as a fixed-amplitude
    drive, the physical picture for linearly driven modes whose saturation lies
    outside this network. amp_cap (in units of max|y0|) stops the integration
    when any amplitude exceeds it -- unbounded parent growth otherwise shrinks
    the step size without limit and the run never returns."""

    def f(t, y):
        z = rhs(t, y[:4] + 1j * y[4:], w, gam, k1, k2, d1, d2)
        if pump:
            z[:2] = 0.0
        return np.concatenate([z.real, z.imag])

    y0 = np.asarray(y0, dtype=complex)
    scale = float(np.max(np.abs(y0)))
    events = None
    if amp_cap is not None:
        def blowup(t, y):
            return amp_cap * scale - np.max(np.hypot(y[:4], y[4:]))
        blowup.terminal = True
        events = blowup
    sol = solve_ivp(f, (0.0, t_end), np.concatenate([y0.real, y0.imag]),
                    method="DOP853", t_eval=np.linspace(0.0, t_end, n_out),
                    rtol=rtol, atol=1e-16 * scale, events=events)
    if not sol.success:
        raise RuntimeError(sol.message)
    capped = events is not None and len(sol.t_events[0])
    if capped:
        print(f"  RUNAWAY: amplitudes passed {amp_cap:g} x the drive at "
              f"t = {sol.t_events[0][0] / YR:.3f} yr -- integration stopped")
    return sol.t, sol.y[:4] + 1j * sol.y[4:], bool(capped)


def check_conservation(w, k1, k2) -> None:
    """Undamped, exactly resonant: energy is conserved and the quanta budget
    balances. Nothing else in the suite covers this RHS.

    c belongs to both legs, so the single-triplet Manley-Rowe pairs do not hold
    here. With f1, f2 the fluxes through the two legs,

        dN_P = -f1,  dN_Q = +f1,  dN_c = f1 - f2,  dN_d = +2 f2,

    which leaves N_P + N_Q and N_Q - N_c - N_d/2 invariant.
    """
    gam = np.zeros(4)
    y0 = np.array([2e-5, 1e-5, 3e-6, 1e-6], dtype=complex)
    t, y, _ = integrate(y0, 4e11, w, gam, k1, k2, 0.0, 0.0, n_out=400, rtol=1e-12)
    E = np.abs(y) ** 2
    N = E / np.array(w)[:, None]
    for name, x in (("energy", E.sum(axis=0)),
                    ("N_P + N_Q", N[0] + N[1]),
                    ("N_Q - N_c - N_d/2", N[1] - N[2] - 0.5 * N[3])):
        print(f"  {name:<20s} drift {np.ptp(x) / abs(x[0]):.2e}")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--csv", default="out/four_mode_dsct_M2.0.csv")
    p.add_argument("--out", default="out/m6_network.csv")
    p.add_argument("--q-parent", type=float, default=None,
                   help="parent amplitude; default is 3x the instability threshold")
    p.add_argument("--full", action="store_true",
                   help="let the parents evolve (kappa-mechanism growth included) "
                        "instead of holding them as a fixed pump; stops on runaway")
    a = p.parse_args()

    df = pd.read_csv(a.csv)
    # Rank on the detuned parametric threshold, not on Gamma/gamma_d: the latter
    # drops the [1 + Delta^2/(2 gamma_d)^2] factor, which is ~(89)^2 at the median.
    df["E_c_over_E_th"] = df.q_c**2 / df.E_th_over_E_star
    df = df.sort_values("E_c_over_E_th", ascending=False)
    n_above = int((df.E_c_over_E_th > 1).sum())
    print(f"{len(df)} candidates, {n_above} above threshold at q = {1e-6:g}")
    r = df.iloc[0]

    print(f"\nsystem: ({r.l_a},{r.n_a:+.0f}) + ({r.l_b},{r.n_b:+.0f}) -> "
          f"({r.l_c},{r.n_c:+.0f}) -> 2x({r.l_d},{r.n_d:+.0f})  [{r.comb}]")
    print(f"  f = {r.f_a:.4f}, {r.f_b:.4f}, {r.f_c:.4f}, {r.f_d:.4f} c/d")
    print(f"  kappa  direct {r.kappa_direct:+.3f}   param {r.kappa_param:+.3f}")
    print(f"  Delta  direct {r.delta_direct / CD:+.3e}  param {r.delta_param / CD:+.3e} c/d")
    print(f"  E_c/E_th {r.E_c_over_E_th:.3e} at q = 1e-6 "
          f"(Gamma/gamma_d {r.growth_over_damping:.3e} without detuning)")

    lo, hi = (r.f_a, r.f_b) if r.f_a < r.f_b else (r.f_b, r.f_a)
    g_lo, g_hi = (r.gamma_a, r.gamma_b) if r.f_a < r.f_b else (r.gamma_b, r.gamma_a)
    w = np.array([hi, lo, r.f_c, r.f_d]) * CD
    gam = np.array([g_hi, g_lo, r.gamma_c, r.gamma_d])
    k1, k2 = float(r.kappa_direct), float(r.kappa_param)

    print("\nconservation check (undamped, resonant):")
    check_conservation(w, k1, k2)

    # E_c scales as q^4, so the threshold amplitude follows the quarter power.
    # Driving above it automatically puts Gamma above |Delta_p|/2, which is what
    # makes the integration affordable at all.
    q_thr = 1e-6 * r.E_c_over_E_th ** -0.25
    q0 = a.q_parent if a.q_parent else 3.0 * q_thr
    print(f"\nthreshold parent amplitude {q_thr:.3e}; driving at {q0:.3e}")

    gamma_par_driven = r.gamma_par * (q0 / 1e-6) ** 2
    d1, d2 = float(r.delta_direct), float(r.delta_param)
    y0 = np.array([q0, q0, 1e-3 * q0, 1e-4 * q0], dtype=complex)

    if a.full:
        # Parents grow at |gamma_P| (e-fold ~12 d for the best system), so the
        # amplitudes leave the perturbative regime long before any limit cycle:
        # the physical outcome is a runaway, detected by the cap, not integrated
        # through.
        t_end = 60.0 / max(abs(gam).max(), gamma_par_driven)
        print(f"\nfull system, t_end {t_end / YR:.2f} yr, runaway cap 30x drive")
        t, y, capped = integrate(y0, t_end, w, gam, k1, k2, d1, d2, amp_cap=30.0)
    else:
        # Fixed-amplitude pump: P, Q pinned at q0. Real parents saturate at the
        # observed amplitude by physics outside this network, so pinning them
        # isolates the c <-> d parametric dynamics. The slowest scale is c's
        # off-resonant relaxation 1/|Delta_direct + i gamma_c|.
        t_end = 8.0 / min(np.hypot(d1, gam[2]), gam[3], gamma_par_driven)
        print(f"\npumped system (P, Q frozen), t_end {t_end / YR:.2f} yr")
        t, y, capped = integrate(y0, t_end, w, gam, k1, k2, d1, d2, pump=True)
    E = np.abs(y) ** 2

    print(f"\nintegrated {t[-1] / YR:.3e} yr")
    for i, nm in enumerate("PQcd"):
        print(f"  E_{nm}/E_star  start {E[i, 0]:.3e}  end {E[i, -1]:.3e}  max {E[i].max():.3e}")

    half = len(t) // 2
    if capped:
        print("\n  runaway: no saturated state to analyse (rerun without --full "
              "for the pumped limit cycle)")
    else:
        swing = np.ptp(E[3, half:]) / max(E[3, half:].mean(), 1e-300)
        print(f"\n  daughter energy swing over the second half: {swing:.3f}")
        print("  bounded limit cycle" if swing > 1e-3 and E[3, -1] > E[3, 0]
              else "  fixed point or decay")

        # Slaving is a sub-threshold statement: above E_th the parametric leg
        # drains c and |q_c| sits below the slaved value by construction.
        # Verify it where it should hold, at half the threshold amplitude,
        # where d decays and c relaxes cleanly. The response is s_c mu q_P q_Q
        # with s_c = 2 (distinct parents) -- the bare mu misses the
        # combinatorial factor of the amplitude equation.
        mu = 2.0 * abs(w[2] * k1) / np.hypot(d1, gam[2])
        q_sub = 0.5 * q_thr
        y0s = np.array([q_sub, q_sub, 1e-3 * q_sub, 1e-4 * q_sub], dtype=complex)
        ts, ys, _ = integrate(y0s, max(8.0 / np.hypot(d1, gam[2]), 6.0 / gam[2]),
                              w, gam, k1, k2, d1, d2, pump=True)
        hs = len(ts) // 2
        err = np.abs(np.abs(ys[2, hs:]) - mu * q_sub**2) / (mu * q_sub**2)
        print(f"  q_c vs 2 mu q_P q_Q at 0.5x threshold (sub-threshold, slaved): "
              f"median {np.median(err):.3e}")
        print(f"  E_d there: start {abs(ys[3, 0])**2:.3e} -> end {abs(ys[3, -1])**2:.3e} "
              f"(decays, as it must below threshold)")

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"t_yr": t / YR, **{f"E_{nm}": E[i] for i, nm in enumerate("PQcd")}}).to_csv(
        out, index=False)
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
