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
"""

from __future__ import annotations

import argparse
import pathlib
import re

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
    modes, trip, eth, skipped = {}, [], {}, False
    for line in open(path):
        if m := re.match(r"#\s+(\(\d+,[+-]\d+\)) -> ", line):   # selecting pair: ... E_th E_eq|none
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
            modes[int(f[0])] = dict(gen=int(f[1]), n=int(float(f[2])), l=int(float(f[3])), w=float(f[5]))
        elif len(f) in (4, 5):
            trip.append({int(f[0]), int(f[1]), int(f[2])})
    return modes, trip, eth


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

    # (l, n) and role of each mode from the network file, when summary.csv names one
    modes, eth = {}, {}
    if (a.data / "summary.csv").exists():
        net = pathlib.Path(pd.read_csv(a.data / "summary.csv")["file"].iloc[0])
        if net.exists():
            modes, _, eth = read_network(net)
    lab = {n: f"({modes[int(n)]['l']},{modes[int(n)]['n']:+d})" if int(n) in modes else f"mode {n}" for n in names}
    par = [n for n in names if modes.get(int(n), {}).get("gen") == 0] if modes else names[:a.parents]
    dau = [n for n in names if n not in par]

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
    head = a.title or str(a.data)
    ax.set_title(f"{head}\nparents {', '.join(lab[n] for n in par)};  {len(dau)} daughters", fontsize=10)

    ax.set(yscale="log", xlabel="Time [s]", ylim=(1e-34,1e-1), ylabel=r"Mode Energy  $[E_\star]$")
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"network_run.{ext}")
    plt.close(fig)
    print(f"-> {out}/network_run.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
