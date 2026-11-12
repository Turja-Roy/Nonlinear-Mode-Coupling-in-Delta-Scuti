#!/usr/bin/env python3
"""Two runs of one ladder side by side: mode energies on top, driving and damping power 2|gamma|E
below, same y axes. Saved as steps_<run a>_<run b>.png in the ladder folder.

    python3 scripts/step_compare.py out/networks/dsct_M2.0/ladder_1p/l1_n-1/{d6,d8}
"""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt

from network_plot import draw


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", type=pathlib.Path, nargs=2)
    a = ap.parse_args()
    fig, axes = plt.subplots(2, 2, figsize=(15, 8.5), height_ratios=(3, 2), sharex="col", layout="constrained")
    for j, d in enumerate(a.runs):
        draw(d, f"{d.parents[2].name} {d.parent.parent.name}/{d.parent.name}/{d.name}", axes=axes[:, j])
    for row in axes:
        lo = min(x.get_ylim()[0] for x in row)
        hi = max(x.get_ylim()[1] for x in row)
        for x in row:
            x.set_ylim(lo, hi)
    out = a.runs[0].parent / f"steps_{a.runs[0].name}_{a.runs[1].name}.png"
    fig.savefig(out)
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
