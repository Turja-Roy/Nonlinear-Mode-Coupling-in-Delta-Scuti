#!/usr/bin/env python3
"""Plot network_run's output: one energy curve per mode, log-y vs t.

    ./build/network_run data/fig6_network.data --out out/network_run
    python3 scripts/network_plot.py --data out/network_run

Reads <data>/energies.csv (t, E_<id>...). A handful of modes get their own
colour and a legend entry; more than that get thin grey lines with no
per-mode legend, the same trick mw25_plot.py's Figure 2 uses for its N
daughters -- a legend with 101 entries is not readable anyway.
"""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt
import pandas as pd

PALETTE = ["#0072B2", "#CC3311", "#000000", "#E69F00", "#33BBEE", "#009E73", "#CC79A7", "#56B4E9"]
MANY = len(PALETTE)  # more modes than this -> thin grey lines, no per-mode legend

plt.rcParams.update({
    "figure.dpi": 150, "savefig.bbox": "tight",
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 10,
})


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

    df = pd.read_csv(a.data / "energies.csv")
    t = df["t"].to_numpy()
    names = [c[2:] for c in df.columns if c != "t"]

    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    if len(names) <= MANY:
        for name, col in zip(names, PALETTE):
            ax.plot(t, df[f"E_{name}"], color=col, lw=1.4, label=f"mode {name}")
        ax.legend(fontsize=8, framealpha=0.9)
        ax.set_title(f"{a.data}", fontsize=10)
    else:
        for name in names[a.parents:]:
            ax.plot(t, df[f"E_{name}"], color="0.4", lw=0.4, alpha=0.5)
        ax.set_title(f"{a.data}  ({len(names)} modes)", fontsize=10)
    for name, col in zip(names[:a.parents], PALETTE[:2] if len(names) > MANY else []):
        ax.plot(t, df[f"E_{name}"], color=col, lw=1.6, label=f"parent {name}")
    if a.parents and len(names) > MANY:
        ax.plot([], [], color="0.4", lw=0.8, label="daughters")
        ax.legend(fontsize=8, framealpha=0.9)
    if a.title:
        ax.set_title(a.title, fontsize=10)

    ax.set(yscale="log", xlabel="Time [s]", ylim=(1e-34,1e-1), ylabel=r"Mode Energy  $[E_\star]$")
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"network_run.{ext}")
    plt.close(fig)
    print(f"-> {out}/network_run.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
