#!/usr/bin/env python3
"""MW25 four-mode mixed systems on the wide net.

    a, b parents  ->  c daughter, resonant with omega_a +- omega_b  (direct leg)
    c             ->  d, d granddaughters, resonant with omega_c/2  (parametric leg)

Parent means self-excited (gamma < 0), as in MW23 sec 2. c is damped and born
from a and b, so its own decay products d are granddaughters (MW23 sec 5).

Ranking runs on kappa at m = (0, 0, 0): angular_factors goes through sympy's
wigner_3j, and the full m set of an l ~ 15 triplet is a few hundred symbols,
which would dominate the cost across ~47k triplets. The top of the ranking then
gets the full m treatment.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import pathlib
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from coupling.angular import satisfies_selection_rules
from coupling.kappa import kappa_abc, kappa_all_m
from coupling.modes import load_model
from coupling.observables import threshold_energy
from coupling.triplets import DETUNING_CUT_DIMLESS, RadialTriplet

CD = 7.27220522e-5  # rad/s per cycle/day
Q_PARENT = 1e-6  # MW23 Fig. 6 reference parent amplitude


def _kappa_one(arg):
    key, t, efs = arg
    r = kappa_abc(*(efs[k] for k in t.keys), (0, 0, 0))
    return key, (r.kappa, r.refine)


def _kappa_shared(item):
    return _kappa_one((item[0], item[1], globals()["_EFS"]))


def parse() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--model", default="models/dsct_M2.0")
    p.add_argument("--detail-dir", default="detail_wide")
    p.add_argument("--inlist", default="gyre_ad_wide.in")
    p.add_argument("--nad", nargs="*", default=["summary_nad_wide.h5", "summary_nad_wide_hi.h5"])
    p.add_argument("--out", default="out/four_mode_dsct_M2.0.csv")
    p.add_argument("--top-full-m", type=int, default=200)
    p.add_argument("-j", "--jobs", type=int, default=mp.cpu_count())
    return p.parse_args()


def _cache_key(key: tuple) -> str:
    return ";".join(f"{l},{n}" for l, n in key)


def _load_cache(path: pathlib.Path) -> dict[tuple, tuple[float, int]]:
    """kappa integrals already done by a previous (possibly crashed) run."""
    if not path.exists():
        return {}
    out = {}
    for line in path.read_text().splitlines()[1:]:
        k, kappa, refine = line.split("\t")
        key = tuple(tuple(int(v) for v in p.split(",")) for p in k.split(";"))
        out[key] = (float(kappa), int(refine))
    return out


def main() -> int:
    a = parse()
    t0 = time.time()
    bg, efs = load_model(a.model, detail_dir=a.detail_dir, inlist=a.inlist, nad=tuple(a.nad))
    efs = {k: e for k, e in efs.items() if np.isfinite(e.gamma)}
    cut = DETUNING_CUT_DIMLESS * bg.omega_dyn

    keys = sorted(efs, key=lambda k: efs[k].omega)
    w = np.array([efs[k].omega for k in keys])
    g = np.array([efs[k].gamma for k in keys])
    L = np.array([k[0] for k in keys])
    driven = np.where(g < 0)[0]
    by_l = {l: np.sort(np.where(L == l)[0]) for l in sorted(set(L.tolist()))}
    for l in by_l:
        by_l[l] = by_l[l][np.argsort(w[by_l[l]])]
    print(f"{len(keys)} modes with gamma, {len(driven)} driven, cut {cut / CD:.3f} c/d")

    direct = []
    for i, ia in enumerate(driven):
        for ib in driven[i:]:
            for wc, sgn in ((w[ia] + w[ib], +1), (abs(w[ia] - w[ib]), -1)):
                if wc <= 0.05 * CD:
                    continue
                for lc in by_l:
                    if not satisfies_selection_rules(int(L[ia]), int(L[ib]), lc):
                        continue
                    idx = by_l[lc]
                    hit = idx[(np.abs(w[idx] - wc) < cut) & (g[idx] > 0)]
                    direct += [(int(ia), int(ib), int(k), sgn) for k in hit]

    quads = []
    for ia, ib, ic, sgn in direct:
        lc = int(L[ic])
        if lc % 2:
            continue
        for ld in by_l:
            if 2 * ld < lc:
                continue
            idx = by_l[ld]
            hit = idx[(np.abs(2 * w[idx] - w[ic]) < cut) & (g[idx] > 0)]
            quads += [(ia, ib, ic, int(k), sgn) for k in hit]
    print(f"{len(direct)} direct triples, {len(quads)} four-mode candidates "
          f"[{time.time() - t0:.0f} s]")
    if not quads:
        print("no candidates")
        return 1

    def triplet(sum_slot: int, pair: tuple[int, int]) -> RadialTriplet:
        return RadialTriplet(
            sum_mode=keys[sum_slot],
            pair=(keys[pair[0]], keys[pair[1]]),
            omega=(-float(w[sum_slot]), float(w[pair[0]]), float(w[pair[1]])),
            delta=float(w[pair[0]] + w[pair[1]] - w[sum_slot]),
        )

    # The sign assignment does not enter kappa (only omega^2 does), so one
    # kappa per unordered (l, n) triple serves every candidate that uses it.
    legs: dict[tuple, RadialTriplet] = {}
    for ia, ib, ic, sgn in direct:
        sum_slot, pair = (ic, (ia, ib)) if sgn > 0 else ((ia, (ib, ic)) if w[ia] > w[ib]
                                                         else (ib, (ia, ic)))
        legs.setdefault(tuple(sorted((keys[ia], keys[ib], keys[ic]))), triplet(sum_slot, pair))
    for _, _, ic, idd, _ in quads:
        legs.setdefault(tuple(sorted((keys[ic], keys[idd], keys[idd]))), triplet(ic, (idd, idd)))
    print(f"{len(legs)} distinct radial triplets to integrate")

    t1 = time.time()
    cache_path = pathlib.Path(a.out).with_suffix(".kappa_cache.tsv")
    kap = _load_cache(cache_path)
    items = [(k, t) for k, t in legs.items() if k not in kap]
    if kap:
        print(f"resuming: {len(kap)} triplets cached, {len(items)} to go")

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "a") as fh:
        if fh.tell() == 0:
            fh.write("key\tkappa\trefine\n")

        def _record(key, val):
            kap[key] = val
            fh.write(f"{_cache_key(key)}\t{val[0]!r}\t{val[1]}\n")
            fh.flush()

        if a.jobs == 1:
            for k, t in items:
                _record(*_kappa_one((k, t, efs)))
        else:
            # fork, so the eigenfunctions are shared copy-on-write rather than
            # pickled. Small chunks keep the checkpoint granularity fine: a
            # crash loses at most jobs * chunksize integrals.
            globals()["_EFS"] = efs
            with mp.get_context("fork").Pool(a.jobs) as pool:
                for k, v in pool.imap_unordered(_kappa_shared, items, chunksize=8):
                    _record(k, v)
    print(f"kappa done, {len(kap)} triplets on {a.jobs} core(s) [{time.time() - t1:.0f} s]")

    rows = []
    for ia, ib, ic, idd, sgn in quads:
        kd, rd = kap[tuple(sorted((keys[ia], keys[ib], keys[ic])))]
        kp, rp = kap[tuple(sorted((keys[ic], keys[idd], keys[idd])))]
        wa, wb, wc, wd = w[ia], w[ib], w[ic], w[idd]
        # The sum slot carries the negative sign: the sum branch puts c there,
        # the difference branch the higher-frequency of a, b.
        delta_d = (wa + wb - wc) if sgn > 0 else (min(wa, wb) + wc - max(wa, wb))
        delta_p = 2 * wd - wc
        # q_c driven off resonance by the two parents, then Gamma of the c -> d,d
        # parametric leg at that q_c. s_c = 2 for distinct parents, 1 for a
        # self-coupled pair (a = b, sum branch only) -- the m6 integrator
        # measured the 2 directly against the bare mu.
        s_c = 2.0 if ia != ib else 1.0
        q_c = s_c * abs(wc * kd) / np.hypot(delta_d, g[ic]) * Q_PARENT**2
        gamma_par = 2 * abs(kp) * wd * q_c
        rows.append(dict(
            l_a=keys[ia][0], n_a=keys[ia][1], l_b=keys[ib][0], n_b=keys[ib][1],
            l_c=keys[ic][0], n_c=keys[ic][1], l_d=keys[idd][0], n_d=keys[idd][1],
            comb="sum" if sgn > 0 else "diff",
            f_a=wa / CD, f_b=wb / CD, f_c=wc / CD, f_d=wd / CD,
            gamma_a=g[ia], gamma_b=g[ib], gamma_c=g[ic], gamma_d=g[idd],
            delta_direct=delta_d, delta_param=delta_p,
            kappa_direct=kd, kappa_param=kp, refine_direct=rd, refine_param=rp,
            q_c=q_c, gamma_par=gamma_par, growth_over_damping=gamma_par / g[idd],
            E_c_over_E_star=q_c**2,
            E_th_over_E_star=(eth := threshold_energy(kp, wd, wd, g[idd], g[idd], delta_p, 1.0)),
            E_c_over_E_th=q_c**2 / eth if eth else np.inf,
        ))

    df = pd.DataFrame(rows).sort_values("E_c_over_E_th", ascending=False)
    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\n{len(df)} rows -> {out}")
    print(df.head(15).to_string(index=False,
          columns=["l_a", "n_a", "l_b", "n_b", "l_c", "n_c", "l_d", "n_d", "comb",
                   "f_c", "f_d", "kappa_direct", "kappa_param", "E_c_over_E_th"]))

    # Gamma > gamma_d is not the criterion: the parametric threshold carries the
    # detuning factor [1 + Delta^2/(gamma_b+gamma_c)^2], and |Delta_p|/2gamma_d
    # has median 89 here, so ranking on Gamma/gamma_d overstates the instability
    # by that factor squared.
    print(f"\n{int((df.E_c_over_E_th > 1).sum())} candidates above threshold "
          f"(E_c > E_th) at q_parent = {Q_PARENT:g}; "
          f"{int((df.growth_over_damping > 1).sum())} if the detuning factor is dropped")

    print(f"\nfull m for the top {a.top_full_m}:")
    seen, full = set(), []
    for r in df.itertuples():
        key = ((r.l_c, r.n_c), (r.l_d, r.n_d))
        if key in seen:
            continue
        seen.add(key)
        t = legs[tuple(sorted((key[0], key[1], key[1])))]
        ks, refine = kappa_all_m(t, efs)
        best = max(ks.items(), key=lambda kv: abs(kv[1]))
        full.append(dict(l_c=r.l_c, n_c=r.n_c, l_d=r.l_d, n_d=r.n_d,
                         kappa_m0=r.kappa_param, kappa_max=best[1], m_max=str(best[0]),
                         n_m=len(ks), refine=refine))
        if len(full) >= a.top_full_m:
            break
    fm = pd.DataFrame(full)
    fm.to_csv(out.with_name(out.stem + "_fullm.csv"), index=False)
    print(fm.head(15).to_string(index=False))
    print(f"\ntotal {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
