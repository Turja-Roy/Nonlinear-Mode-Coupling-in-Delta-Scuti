#!/usr/bin/env python3
"""Direct vs parametric census of the resonant triplets.

    python3 scripts/channels.py --tag dsct_M2.0

Reads observables_<tag>.csv, writes channels_<tag>.csv: one row per coupling
channel with counts at both the radial-triplet and the (triplet, m) level, and
the frequency band of the mode each channel drives.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from coupling.observables import CHANNELS, classify_frame, target_index

CD = 7.27220522e-5  # rad/s per cycle/day
KEYS = ["l_a", "n_a", "l_b", "n_b", "l_c", "n_c"]


def census(df: pd.DataFrame) -> pd.DataFrame:
    df = df.assign(channel=classify_frame(df))
    idx = target_index(df)
    w = df[["omega_a", "omega_b", "omega_c"]].abs().to_numpy() / CD
    rows = np.arange(len(df))

    # Direct channels drive one damped mode; the parametric channel feeds the
    # two damped daughters, so quote their band instead.
    f_lo = f_hi = np.where(idx >= 0, w[rows, np.maximum(idx, 0)], np.nan)
    param = df.channel.to_numpy() == "parametric"
    if param.any():
        f_lo, f_hi = f_lo.copy(), f_hi.copy()
        f_lo[param] = w[param, 1:].min(axis=1)
        f_hi[param] = w[param, 1:].max(axis=1)
    df = df.assign(f_lo=f_lo, f_hi=f_hi,
                   frac_detuning=(df.delta / df.omega_a).abs())

    radial = df.drop_duplicates(KEYS)
    out = []
    for name in CHANNELS:
        m, r = df[df.channel == name], radial[radial.channel == name]
        if not len(m):
            continue
        out.append({
            "channel": name,
            "n_radial": len(r),
            "n_m_rows": len(m),
            "mu_max_median": m.mu_max.median(),
            "mu_max_max": m.mu_max.max(),
            "frac_detuning_median": r.frac_detuning.median(),
            "f_target_q10": np.nanpercentile(m.f_lo, 10) if m.f_lo.notna().any() else np.nan,
            "f_target_q50": np.nanpercentile(m[["f_lo", "f_hi"]], 50) if m.f_lo.notna().any() else np.nan,
            "f_target_q90": np.nanpercentile(m.f_hi, 90) if m.f_hi.notna().any() else np.nan,
        })
    return pd.DataFrame(out).sort_values("n_radial", ascending=False, ignore_index=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="dsct_M2.0")
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("out"))
    args = ap.parse_args()

    src = args.out / f"observables_{args.tag}.csv"
    df = pd.read_csv(src)
    tab = census(df)

    n_radial, n_rows = len(df.drop_duplicates(KEYS)), len(df)
    direct = tab[tab.channel.str.startswith("direct")]
    param = tab[tab.channel == "parametric"]
    print(f"{src}: {n_radial} radial triplets, {n_rows} with m\n")
    with pd.option_context("display.width", 160, "display.max_columns", None):
        print(tab.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(f"\ndirect     {int(direct.n_radial.sum()):6d} radial "
          f"{int(direct.n_m_rows.sum()):7d} with m")
    print(f"parametric {int(param.n_radial.sum()):6d} radial "
          f"{int(param.n_m_rows.sum()):7d} with m")

    dst = args.out / f"channels_{args.tag}.csv"
    tab.to_csv(dst, index=False)
    print(f"-> {dst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
