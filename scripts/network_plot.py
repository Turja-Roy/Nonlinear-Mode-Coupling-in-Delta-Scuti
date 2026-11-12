#!/usr/bin/env python3
"""Plot network_run's output: one energy curve per mode, log-y vs t.

    ./build/network_run data/fig6_network.data --out out/network_run
    python3 scripts/network_plot.py --data out/network_run

Reads <data>/energies.csv (t, E_<id>...). A handful of modes get their own
colour and a legend entry; more than that get thin grey lines with no
per-mode legend, the same trick mw25_plot.py's Figure 2 uses for its N
daughters -- a legend with 101 entries is not readable anyway.

Dashed lines: each parent's lowest E_th over the selecting pairs in the
network file's header, and that pair's E_eq when it has one, in the parent's
colour.

Bottom panel (needs the network file, for gamma): P(t) = sum 2|gamma_p| E_p,
the driving, D(t) = sum 2 gamma_d E_d, the damping, and each mode's own term;
legend: each damped mode's share of the <D> over the second half.
"""

from __future__ import annotations

import argparse
import pathlib
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PALETTE = ["#0072B2", "#CC3311", "#000000", "#E69F00", "#33BBEE", "#009E73", "#CC79A7", "#56B4E9"]
MANY = len(PALETTE)  # more modes than this -> thin grey lines, no per-mode legend
TOTAL = "#882255"  # P and D

plt.rcParams.update({
    "figure.dpi": 150, "savefig.bbox": "tight",
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 10,
})


def read_network(path: pathlib.Path):
    modes, trip, eth, skipped = {}, [], {}, False
    for line in open(path):
        if m := re.match(r"#\s+(?:pair\s+)?(\(\d+,[+-]\d+\)) -> ", line):   # selecting pair: ... E_th E_eq|none
            f = line.split()
            if m[1] not in eth or float(f[-2]) < eth[m[1]][0]:
                eth[m[1]] = (float(f[-2]), None if f[-1] == "none" else float(f[-1]))
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


def power(bx, df, t, names, par, dau, modes, lab) -> None:
    term = {n: 2.0 * abs(modes[int(n)]["g"]) * df[f"E_{n}"].to_numpy() for n in names}
    driven = [n for n in names if modes[int(n)]["g"] < 0.0]
    damped = [n for n in names if n not in driven]
    P = sum(term[n] for n in driven)
    D = sum(term[n] for n in damped)
    half = slice(len(t) // 2, None)
    share = {n: term[n][half].mean() / D[half].mean() for n in damped}
    many = len(names) > MANY
    col = dict(zip(par + dau, PALETTE)) if not many else dict(zip(par, PALETTE))
    if len(driven) > 1:
        for n in driven:
            bx.plot(t, term[n], color=col.get(n, "0.4"), lw=1.0, label=r"$2|\gamma|E$ " + lab[n])
    for n in sorted(damped, key=share.get, reverse=True):
        if n in col:
            bx.plot(t, term[n], color=col[n], lw=1.0, label=rf"$2\gamma E$ {lab[n]}: {100 * share[n]:.0f}% of $\langle D \rangle$")
        else:
            bx.plot(t, term[n], color="0.4", lw=0.4, alpha=0.5)
    if many:
        top = sorted(damped, key=share.get, reverse=True)[:3]
        bx.plot([], [], color="0.4", lw=0.8, label=f"{len(damped)} damped: "
                + ", ".join(f"{lab[n]} {100 * share[n]:.0f}%" for n in top) + " of $\langle D \rangle$")
    bx.plot(t, P, color=TOTAL, lw=2.2, label=r"$P=\sum 2|\gamma_p|E_p$")
    bx.plot(t, D, color=TOTAL, lw=2.2, ls=(0, (4, 2)), label=r"$D=\sum 2\gamma_d E_d$")
    top = max(np.nanmax(P), np.nanmax(D))
    bx.set(yscale="log", ylim=(1e-12 * top, 3.0 * top), xlabel="Time [s]",
           ylabel="Driving / damping power\n" r"$2|\gamma|E$  $[E_\star\,{\rm s}^{-1}]$")
    bx.set_title(f"second half: $\langle D \rangle / \langle P \rangle$ = {D[half].mean() / P[half].mean():.4f}",
                 fontsize=9)
    bx.legend(fontsize=7, framealpha=0.9, ncol=2)


def draw(data: pathlib.Path, title: str | None = None, axes=None, n_par: int = 0):
    """Energies (and, with the network file, 2|gamma|E) of one run, on axes (ax, bx) or a new figure."""
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
        ax, bx = axes
    elif rates:
        fig, (ax, bx) = plt.subplots(2, 1, figsize=(8, 8.5), sharex=True, height_ratios=(3, 2),
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
        power(bx, df, t, names, par, dau, modes, lab)
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
