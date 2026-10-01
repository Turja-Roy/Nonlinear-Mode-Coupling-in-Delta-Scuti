#!/usr/bin/env python3
"""Plot the MW25 figures from the C++ integrator's output.

    ./build/mw25 --out out/mw25_cpp
    python3 scripts/mw25_plot.py --data out/mw25_cpp --figs out/mw25_cpp

apps/mw25.cpp writes one CSV per panel (t plus one energy column per mode) and
a lines.csv of the horizontal thresholds and equilibria. This only draws them.

The physics, the parameter values and the sign-convention translation all live
in apps/mw25.cpp.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Same hue-separated palette as make_figures.py: no two data colours share a
# hue and no pair sits on the red-green axis.
BLUE, YELLOW, RED, CYAN, BLACK = ("#0072B2", "#E69F00", "#CC3311", "#33BBEE", "#000000")
GREY = "0.55"

plt.rcParams.update({
    "figure.dpi": 150, "savefig.bbox": "tight",
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 10,
    "lines.linewidth": 1.4,
})

# mode name -> (colour, label). MW25 draw the parents red and blue, the direct
# daughter black and the parametric daughters grey and purple; purple becomes
# yellow here so the set survives colour blindness.
STYLE = {
    "a":  (BLUE,   r"Parent $a$"),        "b":  (RED,    r"Parent $b$"),
    "c":  (BLACK,  r"Daughter $c$"),      "d":  (GREY,   r"Parametric daughter $d$"),
    "d1": (GREY,   r"Parametric daughter $d_1$"),
    "d2": (YELLOW, r"Parametric daughter $d_2$"),
}
DATA = FIGS = None
LINES: pd.DataFrame | None = None


def load(name):
    df = pd.read_csv(DATA / f"{name}.csv")
    return df["t"].to_numpy(), {c[2:]: df[c].to_numpy() for c in df.columns if c != "t"}


def lines_for(panel):
    if LINES is None:
        return {}
    m = LINES[LINES.panel.astype(str) == panel]
    return dict(zip(m.label, m.value))


def draw(ax, t, E, order=None, lw=1.4):
    for name in (order or E):
        col, lbl = STYLE.get(name, (GREY, name))
        ax.plot(t, E[name], color=col, lw=lw, label=lbl)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIGS / f"{name}.{ext}")
    plt.close(fig)
    print(f"  {name}")


def fig1():
    fig, axes = plt.subplots(2, 1, figsize=(6.4, 7.2), sharex=True, layout="constrained")
    t, E = load("fig1a")
    axes[0].plot(t, E["a"], color=BLUE, label=r"Parent $a = b$")
    axes[0].plot(t, E["c"], color=BLACK, label=r"Daughter $c$")
    # for lbl, v in lines_for("1a").items():          # "E_a_eq" -> $E_{a,\rm eq}$
    #     slot = lbl.split("_")[1]
    #     axes[0].axhline(v, color="0.3", ls=":", lw=1.1,
    #                     label=rf"$E_{{{slot},\rm eq}}$ (Eq. 6)")
    # measured plateau (median of last 10% of t); analytic E_a_eq is negative
    # in the self-coupled case so the lines.csv values land off-canvas.
    tail = slice(-max(1, len(t) // 10), None)
    Ec_m = float(np.median(E["c"][tail]))
    Ea_m = float(np.median(E["a"][tail]))
    for v, lbl in ((Ec_m, r"$E_{c,\rm meas}$"), (Ea_m, r"$E_{a,\rm meas}$")):
        if v > 0:
            axes[0].axhline(v, color="0.3", ls=":", lw=1.1, label=lbl)
    y_top = max(1e-3, 10.0 * max(Ec_m, Ea_m, 0.0))
    axes[0].set(xlim=(0, 1200),
                title="(a) self-coupled parent: stable")

    t, E = load("fig1b")
    draw(axes[1], t, E, order=("a", "b", "c"))
    axes[1].set(ylim=(1e-16, 1e1), title="(b) distinct parents: unstable")
    for ax in axes:
        ax.set(yscale="log", ylabel=r"Mode Energy  $[E_\star]$")
        ax.legend(fontsize=8, framealpha=0.9, loc="lower right")
    axes[1].set_xlabel("Time")
    fig.suptitle("MW25 Figure 1 -- direct coupling alone", fontsize=11)
    save(fig, "mw25_fig1")


def fig2():
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.2), sharey=True, layout="constrained")
    for ax, N in zip(axes, (2, 10, 50)):
        t, E = load(f"fig2_N{N}")
        for j in range(N):
            ax.plot(t, E[f"d{j}"], color=BLACK, lw=0.5, alpha=0.6,
                    label=f"{N} daughters" if j == 0 else None)
        ax.plot(t, E["a"], color=BLUE, label=r"Parent $a$")
        ax.plot(t, E["b"], color=RED, label=r"Parent $b$")
        ax.set(yscale="log", ylim=(1e-22, 1e-2), xlabel="Time", title=f"$N = {N}$")
        ax.legend(fontsize=8, framealpha=0.9, loc="lower right")
    axes[0].set_ylabel(r"Mode Energy  $[E_\star]$")
    fig.suptitle("MW25 Figure 2 -- two parents, $N$ damped daughters, still unstable",
                 fontsize=11)
    save(fig, "mw25_fig2")


def _mixed(ax, name, panel, window, ylim):
    t, E = load(name)
    draw(ax, t, E, order=("a", "b", "c", "d"))
    for lbl, v in lines_for(panel).items():
        ax.axhline(v, color="0.3", ls="--" if lbl == "E_a_th" else ":", lw=1.1,
                   label="parametric threshold" if lbl == "E_a_th" else r"$E_{a,\rm eq}$ (Eq. A7)")
    ax.set(yscale="log", xlim=window, ylim=ylim, xlabel="Time")
    ax.legend(fontsize=7, framealpha=0.9, loc="lower left", ncols=2)


def fig3():
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
    _mixed(ax, "fig3", "3", (1500, 2000), (1e-10, 1e-1))
    ax.set_ylabel(r"Mode Energy  $[E_\star]$")
    ax.set_title("MW25 Figure 3 -- mixed direct + parametric coupling: limit cycle")
    save(fig, "mw25_fig3")


def fig4():
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6), layout="constrained")
    _mixed(axes[0], "fig4a", "4a", (0, 5000), (1e-21, 1e-1))
    axes[0].set_title(r"(a) $\kappa_{\rm param} = 100\,\kappa_{\rm direct}$")
    _mixed(axes[1], "fig4b", "4b", (0, 500000), (1e-21, 1e-1))
    axes[1].set_title(r"(b) all linear rates $\times\, 0.01$")
    axes[0].set_ylabel(r"Mode Energy  $[E_\star]$")
    fig.suptitle("MW25 Figure 4 -- limit-cycle behaviour vs the parameters", fontsize=11)
    save(fig, "mw25_fig4")


def fig5():
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.4), layout="constrained")
    windows = {  # hard-coded limit-cycle windows per panel
        "fig5a": (1500, 2000),
        "fig5b": (35500, 40000),
    }
    for ax, name, panel, ttl in ((axes[0], "fig5a", "5a", "(a) as Figure 3"),
                                 (axes[1], "fig5b", "5b",
                                  r"(b) parent driving $\times\, 0.3$")):
        t, E = load(name)
        m = lines_for(panel).get("mu", np.nan)
        ax.plot(t, E["c"], color=BLACK, label=r"$E_c = |q_c|^2$")
        ax.plot(t, m ** 2 * E["a"] * E["b"], color=GREY, lw=1.8, alpha=0.9,
                label=r"$\mu^2 E_a E_b$  (Eq. 8)")
        # per-panel auto-scale y, hard-coded x window
        x0, x1 = windows[name]
        mask = (t >= x0) & (t <= x1)
        Ec = E["c"][mask]; M2 = (m ** 2 * E["a"] * E["b"])[mask]
        y_lo = max(min(Ec[Ec > 0].min(), M2[M2 > 0].min()) / 5.0, 1e-30)
        y_hi = max(Ec.max(), M2.max()) * 5.0
        ax.set(yscale="log", xlim=(x0, x1), ylim=(y_lo, y_hi),
               xlabel="Time", title=f"{ttl},  $\\mu = {m:.0f}$")
        ax.legend(fontsize=8, framealpha=0.9, loc="lower left")
    axes[0].set_ylabel(r"Energy  $[E_\star]$")
    fig.suptitle(r"MW25 Figure 5 -- $q_c = \mu\, q_a q_b$ holds under mixed coupling",
                 fontsize=11)
    save(fig, "mw25_fig5")


def fig6():
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8), layout="constrained")
    t, E = load("fig6a")
    draw(axes[0], t, E, order=("a", "b", "c"))
    axes[0].set_title("(a) direct triplet only: unstable")

    t, E = load("fig6b")
    draw(axes[1], t, E, order=("a", "b", "c", "d1", "d2"), lw=1.1)
    for lbl, v in lines_for("6b").items():
        axes[1].axhline(v, color="0.3", ls="--" if lbl.endswith("a") else ":", lw=1.1,
                        label=f"threshold, parent ${lbl[-1]}$")
    axes[1].set_title("(b) five-mode mixed network")
    for ax in axes:
        ax.set(yscale="log", ylim=(1e-40, 1e1), xlabel="Time  [days]",
               ylabel=r"Mode Energy  $[E_\star]$")
        ax.legend(fontsize=7, framealpha=0.9, loc="lower right")
    fig.suptitle(r"MW25 Figure 6 -- 2.0 $M_\odot$ $\delta$ Sct parameters", fontsize=11)
    save(fig, "mw25_fig6")


ALL = {"1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5, "6": fig6}


def main() -> int:
    global DATA, FIGS, LINES
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=pathlib.Path, default=pathlib.Path("out/mw25_cpp"))
    ap.add_argument("--figs", type=pathlib.Path, default=None)
    ap.add_argument("--only", nargs="*", choices=sorted(ALL), default=sorted(ALL))
    a = ap.parse_args()
    DATA = a.data
    FIGS = a.figs or a.data
    FIGS.mkdir(parents=True, exist_ok=True)
    lp = DATA / "lines.csv"
    LINES = pd.read_csv(lp) if lp.exists() else None

    print(f"MW25 figures from {DATA} -> {FIGS}")
    for k in a.only:
        try:
            ALL[k]()
        except FileNotFoundError as e:
            print(f"  fig{k}: missing {pathlib.Path(e.filename).name} -- run ./build/mw25")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
