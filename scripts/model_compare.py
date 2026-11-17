#!/usr/bin/env python3
"""Two builds of one stellar model, compared: the MESA profile GYRE read, the mode frequencies and linear
rates, and the eigenfunctions.

    python3 scripts/model_compare.py models/dsct_M2.0_T7696 models/dsct_M2.0_T7696_new

Each model's GYRE input is the `file = ` of its gyre/gyre_ad_wide.in. Per GYRE summary file both have,
modes are matched by (l, n_pg); omega = Re(freq), and Im(freq) of a non-adiabatic summary is the linear rate
gamma_rad (up to GYRE's sign and units, the same for both). Eigenfunctions: xi_r of every detail_wide file
both have, largest |xi_new - xi_old| over largest |xi_old|. Prints the largest differences; nothing written.
"""

from __future__ import annotations

import argparse
import filecmp
import pathlib
import re

import h5py
import numpy as np


def gyre_input(model: pathlib.Path) -> pathlib.Path:
    m = re.search(r"file\s*=\s*'([^']+)'", (model / "gyre" / "gyre_ad_wide.in").read_text())
    return (model / "gyre" / m.group(1)).resolve()


def mesa_header(profile: pathlib.Path) -> dict[str, float | str]:
    with open(profile) as f:
        lines = [next(f) for _ in range(3)]

    def num(x):
        try:
            return float(x)
        except ValueError:
            return x.strip('"')
    return dict(zip(lines[1].split(), (num(x) for x in lines[2].split())))


def summary(path: pathlib.Path) -> dict[tuple[int, int], complex]:
    with h5py.File(path) as f:
        fr = f["freq"][()]
        return {(int(l), int(n)): complex(z["re"], z["im"]) for l, n, z in zip(f["l"][()], f["n_pg"][()], fr)}


def xi_r(path: pathlib.Path) -> np.ndarray:
    with h5py.File(path) as f:
        z = f["xi_r"][()]
        return z["re"] + 1j * z["im"] if z.dtype.names else z


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("old", type=pathlib.Path)
    ap.add_argument("new", type=pathlib.Path)
    a = ap.parse_args()

    po, pn = gyre_input(a.old), gyre_input(a.new)
    print(f"GYRE input: {po}\n            {pn}\n  byte-identical: {filecmp.cmp(po, pn, shallow=False)}")
    ho, hn = (mesa_header(p.with_suffix("")) for p in (po, pn))       # profileN.data beside profileN.data.GYRE
    for k in ("version_number", "MESA_SDK_version", "compiler", "build", "date"):
        print(f"  {k:16s} {ho.get(k)}  {hn.get(k)}")
    for k in ("model_number", "num_zones", "initial_mass", "star_age", "Teff", "photosphere_L", "photosphere_r"):
        if k in ho and k in hn:
            print(f"  {k:16s} {ho[k]:.10g}  {hn[k]:.10g}  rel diff {abs(hn[k] / ho[k] - 1) if ho[k] else 0:.2e}")

    for s in sorted(p.name for p in (a.old / "gyre").glob("summary_*.h5")):
        if not (a.new / "gyre" / s).exists():
            print(f"{s}: missing in {a.new}")
            continue
        o, n = summary(a.old / "gyre" / s), summary(a.new / "gyre" / s)
        both = sorted(set(o) & set(n))
        dw = max(abs(n[k].real / o[k].real - 1) for k in both)
        dg = max(abs(n[k].imag - o[k].imag) / max(abs(o[k].imag), 1e-300) for k in both)
        print(f"{s}: {len(both)} modes in both, {len(set(o) ^ set(n))} in one only; largest rel diff "
              f"Re(freq) {dw:.2e}, Im(freq) {dg:.2e}")

    do, dn = a.old / "gyre" / "detail_wide", a.new / "gyre" / "detail_wide"
    names = sorted(p.name for p in do.glob("detail.*.h5") if (dn / p.name).exists())
    worst, which = 0.0, ""
    for name in names:
        xo, xn = xi_r(do / name), xi_r(dn / name)
        d = np.inf if xo.shape != xn.shape else np.max(np.abs(xn - xo)) / np.max(np.abs(xo))
        if d > worst:
            worst, which = d, name
    print(f"detail_wide: {len(names)} eigenfunctions in both; largest |xi_r diff|/max|xi_r| {worst:.2e} ({which})"
          + ("  (inf: different grids)" if np.isinf(worst) else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
