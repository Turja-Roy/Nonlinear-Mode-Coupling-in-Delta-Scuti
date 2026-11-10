#!/usr/bin/env python3
"""dq_dt against dA_dt, one figure per case, plus table.txt.

    python3 scripts/slow_amplitude_plot.py --data out/slow_amplitude

Top: E(t), q (wide, light) and A (dashed). Bottom: |q_x - q_y| / |q_x|, largest
over modes; grey: |w| t eps, the rounding of the phase w t.
"""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PALETTE = ["#0072B2", "#CC3311", "#009E73", "#E69F00", "#33BBEE", "#CC79A7", "#56B4E9", "#000000"]
PAIR_STYLE = ["-", "--", "-.", ":"]

plt.rcParams.update({
    "figure.dpi": 150, "savefig.bbox": "tight",
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 10,
})


def plot_case(data: pathlib.Path, case: str, out: pathlib.Path) -> None:
    df = pd.read_csv(data / f"{case}.csv")
    modes = pd.read_csv(data / f"{case}_modes.csv", dtype={"mode": str})
    pairs = pd.read_csv(data / f"{case}_pairs.csv")
    t = df["t"].to_numpy()
    toy = case[:1] == "T" and case[1:2].isdigit()

    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, height_ratios=(3, 2), layout="constrained")
    for (_, m), col in zip(modes.iterrows(), PALETTE):
        ax.plot(t, df[f"E_q13_{m['mode']}"], color=col, lw=3.5, alpha=0.3)
        ax.plot(t, df[f"E_A10_{m['mode']}"], color=col, lw=1.1, ls="--", label=f"mode {m['mode']}")
        if np.isfinite(m["E_fixed_point"]):
            ax.axhline(m["E_fixed_point"], color=col, lw=0.9, ls=":")
    ax.plot([], [], color="0.5", lw=3.5, alpha=0.3, label="q, rtol 1e-13")
    ax.plot([], [], color="0.5", lw=1.1, ls="--", label="A, rtol 1e-10")
    if modes["E_fixed_point"].notna().any():
        ax.plot([], [], color="0.5", lw=0.9, ls=":", label="fixed point")
    ax.set(yscale="log", ylabel=r"Mode energy  $[E_\star]$")
    ax.legend(fontsize=8, framealpha=0.9, ncol=2)
    ax.set_title(case.replace("_", " "))

    for (_, p), col, ls in zip(pairs.iterrows(), PALETTE, PAIR_STYLE):
        d = df[[f"d_{p['x']}_{p['y']}_{m}" for m in modes["mode"]]].max(axis=1)
        bx.plot(t, d, color=col, lw=1.2, ls=ls, label=f"{p['x']} vs {p['y']}")
    bx.plot(t, np.abs(modes["omega"]).max() * t * np.finfo(float).eps, color="0.6", lw=0.9,
            label=r"$|\omega|_{\max} t\,\epsilon$")
    bx.set(yscale="log", xlabel="t  [1/|ω| units]" if toy else "t  [s]",
           ylabel=r"$|q_x - q_y| / |q_x|$, largest mode")
    bx.legend(fontsize=8, framealpha=0.9, ncol=2)

    for ext in ("png", "pdf"):
        fig.savefig(out / f"{case}.{ext}")
    plt.close(fig)
    print(f"-> {out}/{case}.png")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=pathlib.Path, default=pathlib.Path("out/slow_amplitude"))
    ap.add_argument("--out", type=pathlib.Path, default=None, help="default: --data")
    a = ap.parse_args()
    out = a.out or a.data
    out.mkdir(parents=True, exist_ok=True)

    cases = sorted(p.name[: -len("_runs.csv")] for p in a.data.glob("*_runs.csv"))
    runs, pairs = [], []
    for c in cases:
        plot_case(a.data, c, out)
        runs.append(pd.read_csv(a.data / f"{c}_runs.csv").assign(case=c))
        pairs.append(pd.read_csv(a.data / f"{c}_pairs.csv").assign(case=c))
    runs = pd.concat(runs)[["case", "run", "rtol", "stop", "t_last", "rows_compared", "steps", "seconds"]]
    pairs = pd.concat(pairs).pivot(index="case", columns=["x", "y"], values="d_max")
    pairs.columns = [f"{x} vs {y}" for x, y in pairs.columns]
    with open(out / "table.txt", "w") as f:
        f.write("Runs\n\n" + runs.to_string(index=False, float_format="%.3g") + "\n\n")
        f.write("Largest |q_x - q_y| / |q_x| over modes and output times\n\n"
                + pairs.to_string(float_format="%.2e") + "\n")
    print((out / "table.txt").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
