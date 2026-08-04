#!/usr/bin/env python3
"""Check the GYRE output before anything downstream consumes it.

    python3 check_modes.py
    python3 check_modes.py --inlist gyre_ad_wide.in --ad summary_ad_wide.h5 \
        --nad summary_nad_wide.h5 --detail detail_wide

Verifies the mode inventory is complete against the windows the inlist actually
asked for, that the scan windows did not clip it, that the damping-rate sign
convention is right, and that dlnrho/dlnr can be reconstructed from As, V_2 and
Gamma_1 (U_D is polytrope-only).
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys

import h5py
import numpy as np


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


def _blocks(text: str, name: str) -> list[dict[str, str]]:
    out = []
    for body in re.findall(rf"&{name}\b(.*?)^\s*/", text, re.S | re.M):
        body = re.sub(r"!.*", "", body)
        out.append(dict(re.findall(r"(\w+)\s*=\s*'?([^'\s,]+)'?", body)))
    return out


def read_windows(inlist: pathlib.Path) -> tuple[dict[int, list[tuple[int, int]]], float, float]:
    """Per-l n_pg windows (a list per l -- 'core' and 'gnet' are disjoint for
    l <= 6) plus the overall frequency scan span."""
    text = inlist.read_text()
    wins: dict[int, list[tuple[int, int]]] = {}
    for b in _blocks(text, "mode"):
        l = int(b.get("l", 0))
        wins.setdefault(l, []).append((int(b.get("n_pg_min", -99999)), int(b.get("n_pg_max", 99999))))
    scans = [(float(b["freq_min"]), float(b["freq_max"])) for b in _blocks(text, "scan")]
    return wins, min(s[0] for s in scans), max(s[1] for s in scans)


def check_inventory(
    d: dict[str, np.ndarray], wins: dict[int, list[tuple[int, int]]], f_lo: float, f_hi: float
) -> bool:
    l, n = d["l"], d["n_pg"]
    freq = np.real(d["freq"])
    ok = True
    print(f"{len(l)} modes total, {len(wins)} l values requested\n")
    print(f"{'l':>3} {'count':>6} {'n_pg range':>14} {'freq range (c/d)':>22} {'gaps':>6}")
    print("-" * 56)
    for li in sorted(wins):
        m = l == li
        if not m.any():
            print(f"{li:>3} {0:>6}   -- no modes found in {wins[li]}")
            ok = False
            continue
        ns = np.sort(n[m])
        want: set[int] = set()
        for lo, hi in wins[li]:
            # Only fault gaps inside what was found: the window edges may be
            # legitimately empty (n_pg = 0 exists for l >= 2 only, and the
            # asymptotic n(omega) law that set the edges is ~1% accurate).
            want |= set(range(max(lo, ns.min()), min(hi, ns.max()) + 1))
        if li == 1:
            want.discard(0)
        gaps = sorted(want - set(ns.tolist()))
        ok &= not gaps
        print(
            f"{li:>3} {m.sum():>6} {ns.min():>6d} ..{ns.max():>5d} "
            f"{freq[m].min():>10.3f} ..{freq[m].max():>9.3f} {len(gaps):>6}"
        )
        if gaps:
            print(f"      missing n_pg: {gaps if len(gaps) < 20 else gaps[:20] + ['...']}")

    lo, hi = freq.min(), freq.max()
    if lo < f_lo * 1.02:
        print(f"\nlowest mode {lo:.3f} c/d sits on the scan floor {f_lo} -- widen it")
        ok = False
    if hi > f_hi * 0.98:
        print(f"\nhighest mode {hi:.3f} c/d sits on the scan ceiling {f_hi} -- widen it")
        ok = False
    return ok


def check_damping(d: dict[str, np.ndarray], expect_driven: bool = True) -> bool:
    omega = d["omega"]
    gamma = -np.imag(omega)
    n, freq = d["n_pg"], np.real(d["freq"])
    driven = gamma < 0

    print(f"\n{driven.sum()} driven modes (gamma < 0) of {len(gamma)}")
    if not expect_driven:
        # A gnet-only run (high-l g-modes near omega_c/2) has no p-modes for the
        # kappa mechanism to excite, so the driven-mode sign test cannot run.
        # Everything should be radiatively damped; more than a few percent
        # driven would mean the sign convention flipped.
        ok = bool(driven.mean() < 0.05)
        if not ok:
            print("  no p-modes requested yet many modes are driven -- sign convention?")
        pos = gamma[gamma > 0]
        if len(pos):
            print(f"  damping rates (dimensionless): {pos.min():.3e} .. {pos.max():.3e}")
        return ok
    if not driven.any():
        print("no driven modes -- the kappa mechanism should excite low-order p-modes here.")
        print("either the sign convention is backwards or the nad run did not converge.")
        return False

    dn, df = n[driven], freq[driven]
    print(f"  driven n_pg range {dn.min()} .. {dn.max()},  freq {df.min():.2f} .. {df.max():.2f} c/d")
    ok = bool((dn > 0).mean() > 0.5)
    if not ok:
        print("  most driven modes are g-modes -- sign convention is probably backwards")
    else:
        print("  mostly p-modes, as expected for the kappa mechanism")

    pos = gamma[gamma > 0]
    print(f"  damping rates (dimensionless): {pos.min():.3e} .. {pos.max():.3e}")
    return ok


def check_coverage(ad: dict[str, np.ndarray], nad: dict[str, np.ndarray], tol: float = 0.3) -> bool:
    """Matched on frequency, exactly as coupling/modes.py does it. GYRE's
    nonadiabatic classifier duplicates and skips n_pg labels for high-order
    g-modes, so a join on (l, n_pg) reports misses that are only mislabels."""
    la, na, fa = ad["l"], ad["n_pg"], np.real(ad["freq"])
    ln, fn = nad["l"], np.real(nad["freq"])

    dup = len(ln) - len({*zip(ln.tolist(), nad["n_pg"].tolist())})
    missing: list[tuple[int, int, float]] = []
    for li in sorted(set(la.tolist())):
        b = np.sort(fn[ln == li])
        for f0, n0 in zip(fa[la == li], na[la == li]):
            if not len(b):
                missing.append((li, int(n0), float(f0)))
                continue
            j = int(np.clip(np.searchsorted(b, f0), 1, len(b) - 1))
            j = j if abs(b[j] - f0) < abs(b[j - 1] - f0) else j - 1
            spacing = (b[min(j + 1, len(b) - 1)] - b[max(j - 1, 0)]) / 2 or abs(f0)
            if abs(b[j] - f0) > tol * spacing:
                missing.append((li, int(n0), float(f0)))

    print(f"  gamma covers {len(la) - len(missing)} of {len(la)} adiabatic modes"
          f" ({dup} duplicate n_pg labels in the nad summary, harmless)")
    for li, n0, f0 in missing:
        print(f"    no gamma: l = {li:2d}  n_pg = {n0:5d}  {f0:.4f} c/d")
    if missing:
        print("  pass require_gamma=True to build_mode_list to drop these")
    # The nad root finder misses a handful of modes at the extreme edge of a
    # window (0.8% on pass 1, 0.9% on pass 2); they are dropped, not chased.
    # Only a miss rate past a few percent means the run itself went wrong.
    return len(missing) <= 0.02 * len(la)


def check_structure(detail: pathlib.Path) -> bool:
    files = sorted(detail.glob("detail.*.h5"))
    if not files:
        print(f"\nno detail files in {detail}")
        return False
    d = load(files[len(files) // 2])
    x, As, V_2, G1, rho = d["x"], d["As"], d["V_2"], d["Gamma_1"], d["rho"]

    # As = dlnP/dlnr / Gamma_1 - dlnrho/dlnr and V = -dlnP/dlnr = V_2 x^2
    dlnrho_dlnr = -V_2 * x**2 / G1 - As

    m = (x > 0.02) & (x < 0.98) & (rho > 0)
    num = np.gradient(np.log(rho[m]), np.log(x[m]))
    err = np.abs(dlnrho_dlnr[m] - num) / np.maximum(np.abs(num), 1e-3)

    print(f"\n{len(files)} detail files, {len(x)} grid points in {files[len(files)//2].name}")
    print("dlnrho/dlnr from (As, V_2, Gamma_1) vs numerical")
    print(f"  median rel err {np.median(err):.2e}, 95th pct {np.percentile(err, 95):.2e}")
    ok = bool(np.median(err) < 1e-3)
    if not ok:
        print("  identity does not hold -- check the As / V_2 sign conventions")

    needed = ["x", "xi_r", "xi_h", "lag_rho", "eul_Phi", "deul_Phi", "rho", "P", "Gamma_1",
              "As", "V_2", "c_1", "M_r"]
    absent = [k for k in needed if k not in d]
    if absent:
        print(f"  detail_item_list is missing {absent} -- coupling/ cannot load these")
        ok = False
    return ok


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--inlist", type=pathlib.Path, default=pathlib.Path("gyre_ad.in"))
    p.add_argument("--ad", type=pathlib.Path, default=pathlib.Path("summary_ad.h5"))
    p.add_argument("--nad", type=pathlib.Path, default=pathlib.Path("summary_nad.h5"))
    p.add_argument("--detail", type=pathlib.Path, default=pathlib.Path("detail"))
    a = p.parse_args()

    if not a.ad.exists():
        print(f"{a.ad} missing; run gyre {a.inlist} first", file=sys.stderr)
        return 1

    wins, f_lo, f_hi = read_windows(a.inlist)
    print(f"=== adiabatic ({a.ad}, windows from {a.inlist}) ===")
    ad = load(a.ad)
    ok = check_inventory(ad, wins, f_lo, f_hi)
    ok &= check_structure(a.detail)

    if a.nad.exists():
        print(f"\n=== nonadiabatic ({a.nad}) ===")
        nad = load(a.nad)
        expect_driven = any(hi > 0 for ws in wins.values() for _, hi in ws)
        ok &= check_damping(nad, expect_driven)
        ok &= check_coverage(ad, nad)
    else:
        print(f"\n{a.nad} missing; run the nonadiabatic inlist for damping rates")

    print("\nOK" if ok else "\nchecks failed")
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
