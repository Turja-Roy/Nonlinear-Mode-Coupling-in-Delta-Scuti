#!/usr/bin/env python3
"""Network observables against daughter count N (or the detuning cut, for a
cut sweep at fixed N), one line per parent starting energy f_p.

    python3 scripts/network_convergence.py --runs out/networks/one_parent

Reads every <runs>/*/summary.csv that network_run wrote, and the
`# network_build ...` header of the network file each one names, for N,
cut and closure. Converged means a curve has flattened in N.
"""

from __future__ import annotations

import argparse
import pathlib
import re

import matplotlib.pyplot as plt
import pandas as pd

PALETTE = ["#0072B2", "#CC3311", "#E69F00", "#009E73", "#CC79A7", "#000000"]
PANELS = [("E_parents", r"$\langle E_{\rm parents}\rangle\ [E_\star]$", True),
          ("swing", "parent energy swing", False),
          ("balance", "dissipation / driving", False),
          ("n_active", "active daughters", False)]

plt.rcParams.update({
    "figure.dpi": 150, "savefig.bbox": "tight",
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 10,
})
HEADER = re.compile(r"daughters (\d+)\s+cut (\S+)\s+closure (\S+)")


def load(runs: pathlib.Path) -> pd.DataFrame:
    rows = []
    for p in sorted(runs.glob("*/summary.csv")):
        r = pd.read_csv(p).iloc[0].to_dict()
        head = next(l for l in open(r["file"]) if l.startswith("# network_build"))
        n, cut, closure = HEADER.search(head).groups()
        r.update(N=int(n), cut=float(cut), closure=closure != "off", run=p.parent.name)
        rows.append(r)
    df = pd.DataFrame(rows)
    df["balance"] = df["dissipation"] / df["driving"]
    if "f_parent" not in df:                     # runs from before --eth-parent
        df["f_parent"] = 0.0
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, default=None)
    a = ap.parse_args()
    df = load(a.runs)
    df.sort_values(["f_parent", "closure", "cut", "N"]).to_csv(a.runs / "convergence.csv", index=False)

    # a cut sweep at fixed N plots against the cut instead
    x, other = ("N", "cut") if df["N"].nunique() > 1 else ("cut", "N")
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.5), layout="constrained", sharex=True)
    fps = sorted(df["f_parent"].unique())
    for (fp, closure, o), g in df.groupby(["f_parent", "closure", other]):
        g = g.sort_values(x)
        label = f"f_p {fp:g}, {other} {o:g}" + ("" if closure else ", no closure")
        style = dict(color=PALETTE[fps.index(fp) % len(PALETTE)], lw=2, marker="o", ms=5,
                     ls="-" if closure else "--", label=label)
        for ax, (col, _, _) in zip(axes.flat, PANELS):
            ax.plot(g[x], g[col], **style)
            bad = g[g["runaway"].astype(str) == "True"]
            ax.plot(bad[x], bad[col], "x", color="k", ms=8)
            part = g[g["stop"].isin(["wall", "failed"])]   # cut short: not a converged value
            ax.plot(part[x], part[col], "s", mfc="none", color="k", ms=9)
    for ax, (_, ylabel, logy) in zip(axes.flat, PANELS):
        ax.set(xscale="log" if x == "cut" else "linear", ylabel=ylabel, yscale="log" if logy else "linear")
        ax.set_xticks(sorted(df[x].unique()), [f"{v:g}" for v in sorted(df[x].unique())])
        ax.xaxis.set_minor_locator(plt.NullLocator())
    for ax in axes[1]:
        ax.set_xlabel("daughters N" if x == "N" else r"detuning cut $[\sqrt{GM/R^3}]$")
    axes[0, 0].axhline(1e-12, color="0.3", lw=1, ls=":")
    axes[0, 0].annotate("observed", (1, 1e-12), xycoords=("axes fraction", "data"),
                        ha="right", va="bottom", fontsize=8, color="0.3")
    axes[0, 0].legend(fontsize=8, framealpha=0.9)
    axes[0, 0].set_title(f"{a.runs.name}   (x = runaway, square = stopped by wall time or failure)",
                         fontsize=9, loc="left")

    out = a.out or a.runs
    for ext in ("png", "pdf"):
        fig.savefig(out / f"network_convergence.{ext}")
    print(f"-> {out}/network_convergence.png, {a.runs}/convergence.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
