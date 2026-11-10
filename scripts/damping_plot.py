#!/usr/bin/env python3
"""MW23 Fig. 4: gamma_rad (red x) and gamma_turb (blue x) against omega, one l.
Open circles: |gamma_rad| of driven modes (MW23 leave them out); black rings:
|gamma_rad + gamma_turb| where it is < 0.

    python3 scripts/damping_plot.py --data out/damping_dsct_M2.0_T7696.csv
"""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt
import pandas as pd

RED, BLUE = "#CC3311", "#0072B2"

plt.rcParams.update({
    "figure.dpi": 150, "savefig.bbox": "tight",
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 10,
})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=pathlib.Path, default=pathlib.Path("out/damping_dsct_M2.0_T7696.csv"))
    ap.add_argument("--l", type=int, default=2, help="degree (default 2, as MW23)")
    ap.add_argument("--title", default=None, help="default: the model tag")
    ap.add_argument("--out", type=pathlib.Path, default=None,
                    help="default out/plots/<tag>/fig_damping_l<l>")
    a = ap.parse_args()
    tag = a.data.stem.removeprefix("damping_")
    out = a.out or pathlib.Path("out/plots") / tag / f"fig_damping_l{a.l}"
    out.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(a.data)
    df = df[(df["l"] == a.l) & df["gamma_rad"].notna()]
    damped, driven = df[df["gamma_rad"] > 0], df[df["gamma_rad"] < 0]
    net = df[df["gamma_rad"] + df["gamma_turb"] < 0]

    fig, ax = plt.subplots(figsize=(6, 4.5), layout="constrained")
    ax.scatter(damped["omega"], damped["gamma_rad"], marker="x", s=22, color=RED, lw=1.0,
               label=r"$\gamma_{\rm rad}$")
    ax.scatter(df["omega"], df["gamma_turb"], marker="x", s=22, color=BLUE, lw=1.0,
               label=r"$\gamma_{\rm turb}$")
    ax.scatter(driven["omega"], -driven["gamma_rad"], marker="o", s=26, facecolor="none", edgecolor=RED,
               lw=1.0, label=r"$|\gamma_{\rm rad}|$, driven ($\gamma_{\rm rad}<0$)")
    ax.scatter(net["omega"], -(net["gamma_rad"] + net["gamma_turb"]), marker="o", s=70, facecolor="none",
               edgecolor="k", lw=0.9, label=r"$|\gamma_{\rm rad}+\gamma_{\rm turb}|$, net driven")
    ax.set(xscale="log", yscale="log", xlim=(5e-5, 6e-3), ylim=(1e-14, 1e-4),
           xlabel=r"$\omega$  [rad s$^{-1}$]", ylabel=r"$\gamma$  [s$^{-1}$]")
    ax.set_title(f"{a.title or tag},  $\\ell = {a.l}$:  {len(driven)} driven, {len(net)} net driven, "
                 f"{len(df)} modes")
    ax.legend(fontsize=8, framealpha=0.9, loc="upper left")
    for ext in ("png", "pdf"):
        fig.savefig(out.with_suffix(f".{ext}"))
    plt.close(fig)
    print(f"-> {out}.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
