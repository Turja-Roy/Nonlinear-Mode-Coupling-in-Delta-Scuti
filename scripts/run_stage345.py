#!/usr/bin/env python3
"""Stages 3-5 for one model: triplets -> kappa -> observables.

    python3 scripts/run_stage345.py --model models/dsct_M2.0

Writes three tables under --out, tagged with the model directory name:

    triplets_<tag>.csv     one row per radial triplet
    kappa_<tag>.csv        one row per (triplet, m), with kappa
    observables_<tag>.csv  the same rows plus mu and E_th, ranked by mu_max

The repository had no script for this -- the existing dsct_M2.0 tables were
produced ad hoc -- so the schemas here were reconstructed from those files and
verified against them column by column:

    kappa.mu        == observables.mu_c        (the daughter's mu)
    kappa.abs_kappa == |kappa|
    kappa.max_abs_n == max(|n_a|, |n_b|, |n_c|)
    observables sorted by mu_max descending; triplets and kappa unsorted

kappa_all_m dominates the runtime (~61 s for dsct_M2.0's 71142 m-combinations
at 11.8 ms/triplet), and observables() calls it already, so both tables are
built from one pass rather than two.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from coupling import observables as obs_mod
from coupling import triplets as trip_mod
from coupling.modes import load_model

DETUNING_CUT = 0.15  # in units of sqrt(GM/R^3), following MW23 sec 3.4

KAPPA_COLS = ["l_a", "n_a", "m_a", "l_b", "n_b", "m_b", "l_c", "n_c", "m_c",
              "omega_a", "omega_b", "omega_c", "delta", "kappa", "refine"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", type=pathlib.Path, default="models/dsct_M2.0")
    ap.add_argument("--tag", default=None, help="defaults to the model dir name")
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("out"))
    ap.add_argument("--l-max", type=int, default=3)
    ap.add_argument("--cut", type=float, default=DETUNING_CUT,
                    help="detuning cut in units of sqrt(GM/R^3)")
    ap.add_argument("--force", action="store_true",
                    help="overwrite existing tables for this tag")
    args = ap.parse_args()

    model = pathlib.Path(args.model)
    tag = args.tag or model.name
    args.out.mkdir(parents=True, exist_ok=True)
    paths = {n: args.out / f"{n}_{tag}.csv"
             for n in ("triplets", "kappa", "observables")}
    clash = [p for p in paths.values() if p.exists()]
    if clash and not args.force:
        raise SystemExit("refusing to overwrite: "
                         + ", ".join(str(p) for p in clash) + "  (pass --force)")

    t0 = time.time()
    print(f"loading {model} ...")
    bg, efs = load_model(model)
    cut = args.cut * bg.omega_dyn
    print(f"  {len(efs)} modes, sqrt(GM/R^3) = {bg.omega_dyn:.6e} rad/s, "
          f"detuning cut = {cut:.6e} rad/s")

    print("stage 3: enumerating triplets ...")
    trips = trip_mod.enumerate_triplets(efs, cut, l_max=args.l_max)
    n_m = trip_mod.count_with_m(trips)
    print(f"  {len(trips)} radial triplets, {n_m} with m  [{time.time() - t0:.1f}s]")
    trip_mod.to_frame(trips).to_csv(paths["triplets"], index=False)
    print(f"  -> {paths['triplets']}")

    print("stages 4-5: kappa, mu, thresholds ...")
    df = obs_mod.to_frame(trips, efs, bg)
    print(f"  {len(df)} rows  [{time.time() - t0:.1f}s]")
    df.to_csv(paths["observables"], index=False)
    print(f"  -> {paths['observables']}")

    kap = df[KAPPA_COLS].copy()
    kap["max_abs_n"] = df[["n_a", "n_b", "n_c"]].abs().max(axis=1)
    kap["abs_kappa"] = df["kappa"].abs()
    kap["mu"] = df["mu_c"]
    kap.to_csv(paths["kappa"], index=False)
    print(f"  -> {paths['kappa']}")

    print(f"done in {time.time() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
