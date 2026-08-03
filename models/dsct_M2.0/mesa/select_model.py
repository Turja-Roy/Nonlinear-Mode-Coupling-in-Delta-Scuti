#!/usr/bin/env python3
"""Pick the target profile from LOGS/ and check it.

Selection is on log g, which fixes R and hence E_star and sqrt(GM/R^3).
Teff is reported but not asserted on.

    python3 select_model.py [--logs LOGS] [--logg 3.900] [--by logg|teff]
"""

from __future__ import annotations

import argparse
import math
import pathlib
import sys

G = 6.67430e-8
MSUN = 1.988409870698051e33
RSUN = 6.957e10
LSUN = 3.828e33
DAY = 86400.0

TARGET_LOGG = 3.900
TARGET_TEFF = 7696.0

EXPECT = {"R_Rsun": 2.627, "logg": 3.900, "dyn_freq_cpd": 2.866, "logL": 1.314}
TOL = {"R_Rsun": 0.02, "logg": 0.005, "dyn_freq_cpd": 0.04, "logL": 0.03}


def read_profile_header(path: pathlib.Path) -> dict[str, float]:
    with path.open() as fh:
        fh.readline()
        names = fh.readline().split()
        values = fh.readline().split()
    out: dict[str, float] = {}
    for name, value in zip(names, values):
        try:
            out[name] = float(value.replace("D", "E"))
        except ValueError:
            pass
    return out


def summarize(hdr: dict[str, float]) -> dict[str, float]:
    M = hdr["star_mass"] * MSUN
    R = hdr["photosphere_r"] * RSUN
    L = hdr["photosphere_L"] * LSUN
    return {
        "model": hdr.get("model_number", float("nan")),
        "Teff": hdr["Teff"],
        "R_Rsun": R / RSUN,
        "logL": math.log10(L / LSUN),
        "logg": math.log10(G * M / R**2),
        "dyn_freq_cpd": math.sqrt(G * M / R**3) * DAY / (2 * math.pi),
        "center_h1": hdr.get("center_h1", float("nan")),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--logs", default="LOGS", type=pathlib.Path)
    ap.add_argument("--logg", default=TARGET_LOGG, type=float)
    ap.add_argument("--teff", default=TARGET_TEFF, type=float)
    ap.add_argument("--by", default="logg", choices=("logg", "teff"))
    ap.add_argument("--show", default=8, type=int)
    args = ap.parse_args()

    profiles = sorted(args.logs.glob("profile*.data"))
    if not profiles:
        print(f"no profile*.data in {args.logs}", file=sys.stderr)
        return 1

    rows = []
    for path in profiles:
        try:
            rows.append((path, summarize(read_profile_header(path))))
        except (KeyError, ValueError) as exc:
            print(f"  skipping {path.name}: {exc}", file=sys.stderr)

    if args.by == "logg":
        rows.sort(key=lambda r: abs(r[1]["logg"] - args.logg))
        print(f"matching on log g = {args.logg:.3f}   ({len(rows)} profiles)\n")
    else:
        rows.sort(key=lambda r: abs(r[1]["Teff"] - args.teff))
        print(f"matching on Teff = {args.teff:.0f} K   ({len(rows)} profiles)\n")

    hdr = f"{'profile':>16} {'model':>7} {'Teff':>8} {'R/Rsun':>8} {'logL':>7} {'logg':>7} {'nu_dyn':>8} {'X_c':>7}"
    print(hdr)
    print("-" * len(hdr))
    for path, s in rows[: args.show]:
        print(
            f"{path.name:>16} {int(s['model']):>7d} {s['Teff']:>8.1f} {s['R_Rsun']:>8.3f} "
            f"{s['logL']:>7.3f} {s['logg']:>7.3f} {s['dyn_freq_cpd']:>8.3f} {s['center_h1']:>7.4f}"
        )

    best_path, best = rows[0]
    gyre_path = best_path.with_suffix(".data.GYRE")

    print(f"\nselected: {best_path}")
    print(f"GYRE input: {gyre_path}" + ("" if gyre_path.exists() else "   <-- MISSING"))

    print()
    ok = True
    for key, expected in EXPECT.items():
        got = best[key]
        good = abs(got - expected) <= TOL[key]
        ok &= good
        print(f"  {key:<14} {got:>9.3f}   expect {expected:>9.3f} +/- {TOL[key]:<6.3f}  {'OK' if good else 'OFF'}")

    print(f"\n  Teff {best['Teff']:.1f} K, {best['Teff'] - args.teff:+.1f} K vs MW23")
    print(f"  detuning cut 0.15*nu_dyn = {0.15 * best['dyn_freq_cpd']:.3f} c/d")

    if not ok:
        print("\nchecks failed; do not run GYRE on this model")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
