#!/usr/bin/env python3
"""dq_dt against dA_dt: one figure per case folder written by slow_amplitude, plus table.txt.

    python3 scripts/slow_amplitude_plot.py --data out/q_vs_A

Top: E(t) integrated as dq_dt. Middle: the same modes and start, integrated as dA_dt,
same axes. Bottom: |E_q - E_A| / E_A per mode. Modes labelled (l,n) from the case's
network.data; dotted: the fixed point of a single triplet.
"""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from network_plot import PALETTE, read_network

Q, A = "q13", "A10"


def plot_case(d: pathlib.Path, title: str) -> None:
    df = pd.read_csv(d / "q_vs_A.csv")
    modes = pd.read_csv(d / "modes.csv", dtype={"mode": str})
    rtol = pd.read_csv(d / "runs.csv").set_index("run")["rtol"]
    t = df["t"].to_numpy()
    names = list(modes["mode"])
    toy = not (d / "network.data").exists()
    lab = {n: f"mode {n}" for n in names}
    if not toy:
        net = read_network(d / "network.data")[0]
        lab = {n: ("parent " if net[int(n)]["gen"] == 0 else "") + f"({net[int(n)]['l']},{net[int(n)]['n']:+d})"
               for n in names}

    fig, (ax, bx, cx) = plt.subplots(3, 1, figsize=(8, 10), sharex=True, height_ratios=(3, 3, 1.3),
                                     layout="constrained")
    bx.sharey(ax)
    for x, form, eq in ((ax, Q, "dq_dt"), (bx, A, "dA_dt")):
        for n, col, ef in zip(names, PALETTE, modes["E_fixed_point"]):
            x.plot(t, df[f"E_{form}_{n}"], color=col, lw=1.2, label=lab[n])
            if np.isfinite(ef):
                x.axhline(ef, color=col, lw=0.8, ls=":")
        x.set(yscale="log", ylabel=r"Mode energy  $[E_\star]$")
        x.set_title(f"{eq},  rtol {rtol[form]:.0e}", fontsize=9)
    ax.legend(fontsize=8, framealpha=0.9, ncol=2)
    for n, col in zip(names, PALETTE):
        cx.plot(t, np.abs(df[f"E_{Q}_{n}"] - df[f"E_{A}_{n}"]) / df[f"E_{A}_{n}"], color=col, lw=0.9)
    cx.set(yscale="log", xlabel="t  [1/|ω| units]" if toy else "t  [s]",
           ylabel=r"$|E_q - E_A| \,/\, E_A$")
    fig.suptitle(title, fontsize=10)

    for ext in ("png", "pdf"):
        fig.savefig(d / f"q_vs_A.{ext}")
    plt.close(fig)
    print(f"-> {d}/q_vs_A.png")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=pathlib.Path, default=pathlib.Path("out/q_vs_A"))
    a = ap.parse_args()

    runs, pairs = [], []
    for d in sorted(p.parent for p in a.data.rglob("q_vs_A.csv")):
        c = d.relative_to(a.data).as_posix()
        c = d.name if c == "." else c
        plot_case(d, c)
        runs.append(pd.read_csv(d / "runs.csv").assign(case=c))
        pairs.append(pd.read_csv(d / "pairs.csv").assign(case=c))
    runs = pd.concat(runs)[["case", "run", "rtol", "stop", "t_last", "rows_compared", "steps", "seconds"]]
    pairs = pd.concat(pairs).pivot(index="case", columns=["x", "y"], values="d_max")
    pairs.columns = [f"{x} vs {y}" for x, y in pairs.columns]
    with open(a.data / "table.txt", "w") as f:
        f.write("Runs\n\n" + runs.to_string(index=False, float_format="%.3g") + "\n\n")
        f.write("Largest |q_x - q_y| / |q_x| over modes and output times\n\n"
                + pairs.to_string(float_format="%.2e") + "\n")
    print((a.data / "table.txt").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
