#!/usr/bin/env python3
"""Check the GYRE output before anything downstream consumes it.

    python3 check_modes.py

Verifies the mode inventory is complete, that the scan windows did not clip it,
that the damping-rate sign convention is right, and that dlnrho/dlnr can be
reconstructed from As, V_2 and Gamma_1 (U_D is polytrope-only).
"""

from __future__ import annotations

import pathlib
import sys

import h5py
import numpy as np

AD = pathlib.Path("summary_ad.h5")
NAD = pathlib.Path("summary_nad.h5")
DETAIL = pathlib.Path("detail")

SCAN_MIN, SCAN_MAX = 0.8, 95.0
N_MIN, N_MAX = -20, 20


def load(path: pathlib.Path) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    with h5py.File(path, "r") as f:
        for k, v in f.attrs.items():
            out[k] = v
        for k in f:
            out[k] = f[k][...]
    for k, v in list(out.items()):
        if isinstance(v, np.ndarray) and v.dtype.names == ("re", "im"):
            out[k] = v["re"] + 1j * v["im"]
    return out


def check_inventory(d: dict[str, np.ndarray]) -> bool:
    l, n = d["l"], d["n_pg"]
    freq = np.real(d["freq"])
    ok = True
    print(f"{len(l)} modes total\n")
    print(f"{'l':>3} {'count':>6} {'n_pg range':>14} {'freq range (c/d)':>22} {'gaps':>6}")
    print("-" * 56)
    for li in sorted(set(l)):
        m = l == li
        ns = np.sort(n[m])
        want = set(range(max(N_MIN, ns.min()), min(N_MAX, ns.max()) + 1))
        # l = 1 has no f mode; l >= 2 does, at n_pg = 0
        if li == 1:
            want.discard(0)
        gaps = sorted(want - set(ns.tolist()))
        ok &= not gaps
        print(
            f"{li:>3} {m.sum():>6} {ns.min():>6d} ..{ns.max():>5d} "
            f"{freq[m].min():>10.3f} ..{freq[m].max():>9.3f} {len(gaps):>6}"
        )
        if gaps:
            print(f"      missing n_pg: {gaps}")

    lo, hi = freq.min(), freq.max()
    if lo < SCAN_MIN * 1.02:
        print(f"\nlowest mode {lo:.3f} c/d sits on the scan floor {SCAN_MIN} -- widen it")
        ok = False
    if hi > SCAN_MAX * 0.98:
        print(f"\nhighest mode {hi:.3f} c/d sits on the scan ceiling {SCAN_MAX} -- widen it")
        ok = False

    reached = (n.min() <= N_MIN) and (n.max() >= N_MAX)
    if not reached:
        print(f"\nn_pg spans {n.min()}..{n.max()}, wanted {N_MIN}..{N_MAX}")
    return ok


def check_damping(d: dict[str, np.ndarray]) -> bool:
    omega = d["omega"]
    gamma = -np.imag(omega)
    l, n, freq = d["l"], d["n_pg"], np.real(d["freq"])
    driven = gamma < 0

    print(f"\n{driven.sum()} driven modes (gamma < 0) of {len(gamma)}")
    if not driven.any():
        print("no driven modes -- the kappa mechanism should excite low-order p-modes here.")
        print("either the sign convention is backwards or the nad run did not converge.")
        return False

    dn, dl, df = n[driven], l[driven], freq[driven]
    print(f"  driven n_pg range {dn.min()} .. {dn.max()},  freq {df.min():.2f} .. {df.max():.2f} c/d")
    ok = bool((dn > 0).mean() > 0.5)
    if not ok:
        print("  most driven modes are g-modes -- sign convention is probably backwards")
    else:
        print("  mostly p-modes, as expected for the kappa mechanism")

    pos = gamma[gamma > 0]
    print(f"  damping rates (dimensionless): {pos.min():.3e} .. {pos.max():.3e}")
    return ok


def check_structure() -> bool:
    files = sorted(DETAIL.glob("detail.*.h5"))
    if not files:
        print("\nno detail files found")
        return False
    d = load(files[len(files) // 2])
    x, As, V_2, G1, rho = d["x"], d["As"], d["V_2"], d["Gamma_1"], d["rho"]

    # As = dlnP/dlnr / Gamma_1 - dlnrho/dlnr and V = -dlnP/dlnr = V_2 x^2
    dlnrho_dlnr = -V_2 * x**2 / G1 - As

    m = (x > 0.02) & (x < 0.98) & (rho > 0)
    num = np.gradient(np.log(rho[m]), np.log(x[m]))
    err = np.abs(dlnrho_dlnr[m] - num) / np.maximum(np.abs(num), 1e-3)

    print(f"\ndlnrho/dlnr from (As, V_2, Gamma_1) vs numerical, {files[len(files)//2].name}")
    print(f"  median rel err {np.median(err):.2e}, 95th pct {np.percentile(err, 95):.2e}")
    ok = bool(np.median(err) < 1e-3)
    if not ok:
        print("  identity does not hold -- check the As / V_2 sign conventions")
    return ok


def main() -> int:
    ok = True
    if not AD.exists():
        print(f"{AD} missing; run gyre gyre_ad.in first", file=sys.stderr)
        return 1

    print("=== adiabatic ===")
    ad = load(AD)
    ok &= check_inventory(ad)
    ok &= check_structure()

    if NAD.exists():
        print("\n=== nonadiabatic ===")
        ok &= check_damping(load(NAD))
    else:
        print(f"\n{NAD} missing; run gyre gyre_nad.in for damping rates")

    print("\nOK" if ok else "\nchecks failed")
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
