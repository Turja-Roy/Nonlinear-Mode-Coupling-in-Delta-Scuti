#!/usr/bin/env python3
"""Plot network_run's output: one energy curve per mode, log-y vs t.

    ./build/network_run data/fig6_network.data --out out/network_run
    python3 scripts/network_plot.py --data out/network_run

Reads <data>/energies.csv (t, E_<id>...). A handful of modes get their own
colour and a legend entry; more than that get thin grey lines with no
per-mode legend, the same trick mw25_plot.py's Figure 2 uses for its N
daughters -- a legend with 101 entries is not readable anyway.

Dashed lines: each parent's lowest E_th over the selecting pairs in the
network file's header, and that pair's E_eq when gamma_p + gamma_b + gamma_c > 0,
in the parent's colour.

Middle panel (needs the network file, for gamma): P(t) = sum 2|gamma_p| E_p,
the driving, D(t) = sum 2 gamma_d E_d, the damping, and each mode's own term;
legend: each damped mode's share of <D>_T over the last window.
Bottom panel: <D>_T / <P>_T, the moving time averages over a window
T = WINDOW t_stop, (1/T) * integral over (t - T, t], trapezoid rule on the outputs.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PALETTE = ["#0072B2", "#CC3311", "#000000", "#E69F00", "#33BBEE", "#009E73", "#CC79A7", "#56B4E9"]
MANY = len(PALETTE)  # more modes than this -> thin grey lines, no per-mode legend
TOTAL = "#882255"  # P and D
WINDOW = 0.1  # moving-average window T, fraction of the run's length

plt.rcParams.update({
    "figure.dpi": 150, "savefig.bbox": "tight",
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 10,
})


def read_network(path: pathlib.Path):
    modes, trip, eth, skipped = {}, [], {}, False
    for line in open(path):
        if m := re.match(r"#\s+(?:pair\s+)?(\(\d+,[+-]\d+\)) -> ", line):
            f = line[line.rindex(")") + 1:].split()    # Delta, Delta/omega, kappa, E_th, E_eq|none[, gamma sum]
            stable = f[4] != "none" and (len(f) < 6 or float(f[5]) > 0.0)
            if m[1] not in eth or float(f[3]) < eth[m[1]][0]:
                eth[m[1]] = (float(f[3]), float(f[4]) if stable else None)
        if not line.strip() or line.startswith("#"):
            continue
        if not skipped:
            skipped = True
            continue
        f = line.split()
        if len(f) in (7, 8):
            modes[int(f[0])] = dict(gen=int(f[1]), n=int(float(f[2])), l=int(float(f[3])), w=float(f[5]),
                                    g=float(f[6]))
        elif len(f) in (4, 5):
            trip.append({int(f[0]), int(f[1]), int(f[2])})
    return modes, trip, eth


def window_mean(t, x, w):
    """(t_j, mean of x over the w output intervals ending at t_j), trapezoid rule. A sum of positive terms per
    window: a running integral differenced would cancel once x has dropped by many decades."""
    seg = 0.5 * (x[1:] + x[:-1]) * np.diff(t)
    return t[w:], np.lib.stride_tricks.sliding_window_view(seg, w).sum(axis=1) / (t[w:] - t[:-w])


def powers(df, modes):
    """2|gamma|E of every mode, P (driven modes) and D (damped modes)."""
    term = {n: 2.0 * abs(modes[int(n)]["g"]) * df[f"E_{n}"].to_numpy() for n in (c[2:] for c in df.columns if c != "t")}
    P = sum(v for n, v in term.items() if modes[int(n)]["g"] < 0.0)
    D = sum(v for n, v in term.items() if modes[int(n)]["g"] > 0.0)
    return term, P, D


def window(t) -> int:
    return max(1, round(WINDOW * (len(t) - 1)))


def balance(t, P, D):
    """(t, <D>_T / <P>_T, T), windows of WINDOW of the outputs."""
    w = window(t)
    tt, p = window_mean(t, P, w)
    with np.errstate(divide="ignore", invalid="ignore"):
        return tt, np.where(p > 0.0, window_mean(t, D, w)[1] / p, np.nan), t[-1] - t[-1 - w]


def run_balance(d: pathlib.Path) -> float:
    """<D>_T / <P>_T at the last output of run folder d; NaN without energies.csv or the network file."""
    if not ((d / "energies.csv").exists() and (d / "network.data").exists()):
        return math.nan
    df = pd.read_csv(d / "energies.csv")
    _, P, D = powers(df, read_network(d / "network.data")[0])
    return float(balance(df["t"].to_numpy(), P, D)[1][-1])


def power(bx, cx, df, t, names, par, dau, modes, lab) -> None:
    term, P, D = powers(df, modes)
    driven = [n for n in names if modes[int(n)]["g"] < 0.0]
    damped = [n for n in names if n not in driven]
    tt, R, T = balance(t, P, D)
    w = window(t)
    Dw = window_mean(t, D, w)[1][-1]
    share = {n: window_mean(t, term[n], w)[1][-1] / Dw for n in damped}
    many = len(names) > MANY
    col = dict(zip(par + dau, PALETTE)) if not many else dict(zip(par, PALETTE))
    if len(driven) > 1:
        for n in driven:
            bx.plot(t, term[n], color=col.get(n, "0.4"), lw=1.0, label=r"$2|\gamma|E$ " + lab[n])
    for n in sorted(damped, key=share.get, reverse=True):
        if n in col:
            bx.plot(t, term[n], color=col[n], lw=1.0, label=rf"$2\gamma E$ {lab[n]}: {100 * share[n]:.0f}% of $\langle D \rangle_T$")
        else:
            bx.plot(t, term[n], color="0.4", lw=0.4, alpha=0.5)
    if many:
        top = sorted(damped, key=share.get, reverse=True)[:3]
        bx.plot([], [], color="0.4", lw=0.8, label=f"{len(damped)} damped: "
                + ", ".join(f"{lab[n]} {100 * share[n]:.0f}%" for n in top) + r" of $\langle D \rangle_T$")
    bx.plot(t, P, color=TOTAL, lw=2.2, label=r"$P=\sum 2|\gamma_p|E_p$")
    bx.plot(t, D, color=TOTAL, lw=2.2, ls=(0, (4, 2)), label=r"$D=\sum 2\gamma_d E_d$")
    top = max(np.nanmax(P), np.nanmax(D))
    bx.set(yscale="log", ylim=(1e-12 * top, 3.0 * top),
           ylabel="Driving / damping power\n" r"$2|\gamma|E$  $[E_\star\,{\rm s}^{-1}]$")
    bx.set_title(rf"shares of $\langle D \rangle_T$ over the last window, $T$ = {T:.3g} s", fontsize=9)
    bx.legend(fontsize=7, framealpha=0.9, ncol=2)

    cx.plot(tt, R, color=TOTAL, lw=1.6)
    cx.axhline(1.0, color="0.3", lw=0.8, ls=":")
    lin = np.all((R[np.isfinite(R)] > 0.5) & (R[np.isfinite(R)] < 2.0))
    cx.set(yscale="linear" if lin else "log", xlabel="Time [s]",
           ylabel=r"$\langle D \rangle_T / \langle P \rangle_T$")
    cx.set_title(rf"moving averages over $T$ = {T:.3g} s; at $t_{{\rm stop}}$: {R[-1]:.4g}", fontsize=9)


def draw(data: pathlib.Path, title: str | None = None, axes=None, n_par: int = 0):
    """Energies (and, with the network file, 2|gamma|E and <D>_T/<P>_T) of one run, on axes (ax, bx, cx)
    or a new figure."""
    df = pd.read_csv(data / "energies.csv")
    t = df["t"].to_numpy()
    names = [c[2:] for c in df.columns if c != "t"]

    # the run's own copy of the network file first (network_sweep.sh)
    modes, eth = {}, {}
    net = data / "network.data"
    if not net.exists() and (data / "summary.csv").exists():
        net = pathlib.Path(pd.read_csv(data / "summary.csv")["file"].iloc[0])
    if net.exists():
        modes, _, eth = read_network(net)
    lab = {n: f"({modes[int(n)]['l']},{modes[int(n)]['n']:+d})" if int(n) in modes else f"mode {n}" for n in names}
    par = [n for n in names if modes.get(int(n), {}).get("gen") == 0] if modes else names[:n_par]
    dau = [n for n in names if n not in par]

    rates = bool(modes) and all(int(n) in modes for n in names)
    fig = None
    if axes is not None:
        ax, bx, cx = axes
    elif rates:
        fig, (ax, bx, cx) = plt.subplots(3, 1, figsize=(8, 10), sharex=True, height_ratios=(3, 2, 1.2),
                                         layout="constrained")
    else:
        fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    if len(names) <= MANY:
        for name, col in zip(par + dau, PALETTE):
            ax.plot(t, df[f"E_{name}"], color=col, lw=1.6 if name in par else 1.2,
                    label=("parent " if name in par else "") + lab[name])
    else:
        for name in dau:
            ax.plot(t, df[f"E_{name}"], color="0.4", lw=0.4, alpha=0.5)
        for name, col in zip(par, PALETTE):
            ax.plot(t, df[f"E_{name}"], color=col, lw=1.6, label="parent " + lab[name])
        ax.plot([], [], color="0.4", lw=0.8, label=f"{len(dau)} daughters")
    for name, col in zip(par, PALETTE):
        e_th, e_eq = eth.get(lab[name], (None, None))
        if e_th:
            ax.axhline(e_th, color=col, lw=0.9, ls=(0, (6, 3)), label=r"$E_{\rm th}$ " + lab[name])
        if e_eq:
            ax.axhline(e_eq, color=col, lw=0.9, ls=(0, (2, 2)), label=r"$E_{\rm eq}$ " + lab[name])
    ax.legend(fontsize=8, framealpha=0.9)
    head = title or str(data)
    ax.set_title(f"{head}\nparents {', '.join(lab[n] for n in par)};  {len(dau)} daughters", fontsize=10)

    ax.set(yscale="log", xlabel="" if rates else "Time [s]", ylim=(1e-34,1e-1), ylabel=r"Mode Energy  $[E_\star]$")
    if rates:
        power(bx, cx, df, t, names, par, dau, modes, lab)
    return fig


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=pathlib.Path, default=pathlib.Path("out/network_run"))
    ap.add_argument("--out", type=pathlib.Path, default=None)
    ap.add_argument("--parents", type=int, default=0,
                    help="highlight the first k modes (network_build writes parents first)")
    ap.add_argument("--title", default=None)
    a = ap.parse_args()
    out = a.out or a.data
    fig = draw(a.data, a.title, n_par=a.parents)
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"network_run.{ext}")
    plt.close(fig)
    print(f"-> {out}/network_run.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
