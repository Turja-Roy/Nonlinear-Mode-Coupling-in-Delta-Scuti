#!/usr/bin/env python3
"""Plot network_run's output: one energy curve per mode, log-y vs t.

    ./build/network_run data/fig6_network.data --out out/network_run
    python3 scripts/network_plot.py --data out/network_run

Reads <data>/energies.csv (t, E_<id>...). A handful of modes get their own
colour and a legend entry; more than that get thin grey lines with no
per-mode legend, the same trick mw25_plot.py's Figure 2 uses for its N
daughters -- a legend with 101 entries is not readable anyway.

--copies: each mode is stored twice, +omega and -omega (amplitude.hpp), and
every triplet once per convention. One panel per choice of parent copies
(2 for one parent, 2x2 for two), each with the daughter copies those parent
copies couple to, from <data>/energies_all.csv and the network file that
summary.csv names. A parent copy with no triplet shows no daughters.
"""

from __future__ import annotations

import argparse
import itertools
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


def read_network(path: pathlib.Path):
    modes, trip, skipped = {}, [], False
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        if not skipped:
            skipped = True
            continue
        f = line.split()
        if len(f) in (7, 8):
            modes[int(f[0])] = dict(gen=int(f[1]), n=int(float(f[2])), l=int(float(f[3])), w=float(f[5]))
        elif len(f) in (4, 5):
            trip.append({int(f[0]), int(f[1]), int(f[2])})
    return modes, trip


def copies(a, out: pathlib.Path) -> None:
    df = pd.read_csv(a.data / "energies_all.csv")
    t = df["t"].to_numpy()
    net = pd.read_csv(a.data / "summary.csv")["file"].iloc[0]
    modes, trip = read_network(pathlib.Path(net))
    par = {}                                        # (l, n) -> {+1: id, -1: id}
    for i, m in modes.items():
        if m["gen"] == 0:
            par.setdefault((m["l"], m["n"]), {})[1 if m["w"] > 0 else -1] = i
    keys = sorted(par)
    combos = list(itertools.product((1, -1), repeat=len(keys)))
    nr = 2 if len(keys) == 2 else 1
    fig, axes = plt.subplots(nr, len(combos) // nr, figsize=(12, 4.5 * nr), sharex=True,
                             sharey=True, squeeze=False, layout="constrained")
    for ax, signs in zip(axes.flat, combos):
        sel = [par[k][s] for k, s in zip(keys, signs)]
        dau = sorted({i for tr in trip if tr & set(sel) for i in tr} - {i for p in par.values() for i in p.values()})
        for i in dau:
            ax.plot(t, df[f"E_{i}"], color="0.4", lw=0.5, alpha=0.6)
        for i, k, s, col in zip(sel, keys, signs, PALETTE):
            ax.plot(t, df[f"E_{i}"], color=col, lw=1.6, label=f"({k[0]},{k[1]}) {'+' if s > 0 else '-'}omega, id {i}")
        ax.plot([], [], color="0.4", lw=0.8, label=f"{len(dau)} daughter copies coupled")
        ax.set_title(", ".join(f"({k[0]},{k[1]}) {'+' if s > 0 else '-'}w" for k, s in zip(keys, signs)))
        ax.set(yscale="log", ylim=(1e-34, 1e-1))
        ax.legend(fontsize=7, framealpha=0.9)
    for ax in axes[-1]:
        ax.set_xlabel("Time [s]")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"Mode Energy  $[E_\star]$")
    if a.title:
        fig.suptitle(a.title + "  (omega copies)", fontsize=10)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"network_copies.{ext}")
    plt.close(fig)
    print(f"-> {out}/network_copies.png")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=pathlib.Path, default=pathlib.Path("out/network_run"))
    ap.add_argument("--out", type=pathlib.Path, default=None)
    ap.add_argument("--parents", type=int, default=0,
                    help="highlight the first k modes (network_build writes parents first)")
    ap.add_argument("--title", default=None)
    ap.add_argument("--copies", action="store_true", help="also the +-omega copy panels")
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
    if a.copies:
        copies(a, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
