#!/usr/bin/env python3
"""E_th and E_eq of every candidate triplet, against the parent frequency and against the detuning.

    python3 scripts/energy_frequency.py --csv data/networks/dsct_M2.0_T7696/ladder_1p/candidates.csv \
        --out out/plots/dsct_M2.0_T7696

Reads network_build's candidates.csv. Top row: E_th of the pairs p -> b + c. Bottom row: E_eq of every
triplet; a sum triplet x + y -> c is drawn at f_x and at f_y. Left: against the parent frequency (Weinberg+12
fig. 6 draws the lowest E_th against the orbital period, i.e. the parent frequency). Right: against
|Delta|/omega of the triplet's highest-frequency mode (the parent of a pair, c of a direct triplet; MW23
fig. 5's detuning axis). Filled: gamma sum > 0; open: gamma sum <= 0 (a pair's fixed point unstable, a
harmonic triplet that cannot hold x alone; a sum triplet's fixed point is unstable at any value). Ringed:
each parent's lowest value in the row.
"""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

KIND = {"pair": ("#0072B2", "o", "pair p → b + c"), "harmonic": ("#E69F00", "^", "harmonic x + x → c"),
        "sum": ("#009E73", "s", "sum x + y → c")}

plt.rcParams.update({
    "figure.dpi": 150, "savefig.bbox": "tight",
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 10,
})


def per_parent(df: pd.DataFrame) -> pd.DataFrame:
    """One row per (triplet, parent): a sum triplet once for each of its parents."""
    a = df.assign(p=df["parent"], f=df["f_parent"])
    b = df[df["kind"] == "sum"].assign(p=lambda x: x["parent2"], f=lambda x: x["f_parent2"])
    return pd.concat([a, b], ignore_index=True)


def panel(ax, d: pd.DataFrame, x: str, y: str) -> None:
    d = d[np.isfinite(d[y]) & (d[y] > 0)]
    for k, (col, mk, lab) in KIND.items():
        s = d[d["kind"] == k]
        if s.empty:
            continue
        up = s["gamma_sum"] > 0
        ax.scatter(s.loc[up, x], s.loc[up, y], s=16, marker=mk, color=col, alpha=0.8, lw=0, label=lab)
        ax.scatter(s.loc[~up, x], s.loc[~up, y], s=16, marker=mk, facecolor="none", edgecolor=col, lw=0.8,
                   alpha=0.8, label=lab + r", $\sum\gamma \leq 0$")
    low = d.loc[d.groupby("p")[y].idxmin()]
    ax.scatter(low[x], low[y], s=70, facecolor="none", edgecolor="k", lw=1.0, label="lowest of a parent")
    ax.set_yscale("log")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--title", default=None)
    a = ap.parse_args()
    d = per_parent(pd.read_csv(a.csv))
    d["detuning"] = d["delta_over_omega"].abs()
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharey="row", layout="constrained")
    for i, y in enumerate(("E_th", "E_eq")):
        rows = d[d["kind"] == "pair"] if y == "E_th" else d
        panel(axes[i, 0], rows, "f", y)
        panel(axes[i, 1], rows, "detuning", y)
        axes[i, 0].set_ylabel(rf"$E_{{\rm {y[2:]}}}$  $[E_\star]$")
        axes[i, 1].set_xscale("log")
    for ax in axes[1]:
        ax.set_xlabel(r"parent frequency $f$ [d$^{-1}$]" if ax is axes[1, 0]
                      else r"$|\Delta| / \omega$ (highest-frequency mode of the triplet)")
    axes[0, 0].set_title("pairs: threshold energy", loc="left")
    axes[1, 0].set_title("all triplets: energy at the fixed point (sum triplet: at both parents)", loc="left")
    axes[1, 0].legend(fontsize=8, framealpha=0.9)
    if a.title:
        fig.suptitle(a.title, fontsize=11)
    a.out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(a.out / f"fig_energy_frequency.{ext}")
    print(f"-> {a.out}/fig_energy_frequency.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
