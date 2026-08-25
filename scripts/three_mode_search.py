#!/usr/bin/env python3
"""Three-mode coupling tables for one model: triplets -> kappa -> observables.

    python3 scripts/three_mode_search.py --model models/dsct_M2.0

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

For the wide net the defaults do not scale: point --detail-dir at detail_wide,
restrict the sum mode with --parents-only, and rank on --m000. See
Plans/l-convergence.md.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time
from dataclasses import replace

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from coupling import observables as obs_mod
from coupling import triplets as trip_mod
from coupling.modes import gamma_turb, load_model

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
    ap.add_argument("--detail-dir", default="detail",
                    help="eigenfunction dumps under <model>/gyre; detail_wide "
                         "for the l <= 25 net")
    ap.add_argument("--inlist", default="gyre_ad.in")
    ap.add_argument("--nad", nargs="*", default=["summary_nad.h5"])
    ap.add_argument("--gamma", choices=("rad", "tot"), default="rad",
                    help="tot adds the turbulent rate to the radiative one, "
                         "which also moves the parent/daughter split")
    ap.add_argument("--parents-only", action="store_true",
                    help="restrict the sum mode to the self-excited modes; the "
                         "pair still ranges over the whole net")
    ap.add_argument("--m000", action="store_true",
                    help="kappa at m = (0, 0, 0) only, skipping sympy's "
                         "wigner_3j over every m combination")
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
    bg, efs = load_model(model, detail_dir=args.detail_dir, inlist=args.inlist,
                         nad=tuple(args.nad))
    cut = args.cut * bg.omega_dyn
    print(f"  {len(efs)} modes, sqrt(GM/R^3) = {bg.omega_dyn:.6e} rad/s, "
          f"detuning cut = {cut:.6e} rad/s")

    if args.gamma == "tot":
        n_drv_rad = sum(ef.gamma < 0.0 for ef in efs.values())
        efs = {k: replace(ef, gamma=ef.gamma + gamma_turb(ef)) for k, ef in efs.items()}
        n_drv = sum(ef.gamma < 0.0 for ef in efs.values())
        print(f"  gamma_rad + gamma_turb: driven modes {n_drv_rad} -> {n_drv}")

    sum_keys = None
    if args.parents_only:
        sum_keys = {k for k, ef in efs.items() if ef.gamma < 0.0}
        print(f"  sum mode restricted to {len(sum_keys)} self-excited modes")

    print("enumerating triplets ...")
    trips = trip_mod.enumerate_triplets(efs, cut, l_max=args.l_max, sum_keys=sum_keys)
    print(f"  {len(trips)} radial triplets  [{time.time() - t0:.1f}s]")
    if not args.m000:
        print(f"  {trip_mod.count_with_m(trips)} with m")
    trip_mod.to_frame(trips, with_m=not args.m000).to_csv(paths["triplets"], index=False)
    print(f"  -> {paths['triplets']}")

    print("kappa, mu, thresholds ...")
    df = obs_mod.to_frame(trips, efs, bg, m000=args.m000)
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
