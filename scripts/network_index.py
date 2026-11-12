#!/usr/bin/env python3
"""networks.md for one model: every ladder under <model>/ladder_{1p,2p}/, its parents, how it was
built, its selecting triplets in the order they entered, and one row per run.

    python3 scripts/network_index.py --model out/networks/dsct_M2.0_T7696

Shared daughters and direct triplets are counted from each network's triplet table, not from how it
was built.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import re

import pandas as pd

from network_plot import read_network

MODE = r"\((-?\d+),([+-]?\d+)\)"
NOTATION = """\
Notation. (l,n): degree and radial order of a mode; f: frequency; gamma: linear rate, < 0 driven (parent),
> 0 damped (daughter). N: number of daughters. Pair p -> b + c: parametric triplet; harmonic x + x -> c and
sum x + y -> c: direct triplets. Delta: detuning; kappa: coupling coefficient. E_th: parent energy above
which a pair grows (Weinberg+12 eq. 50); E_eq: parent energy at a triplet's fixed point (MW25 eq. 6), for a
sum triplet E_x + E_y. Energies in E_star = GM^2/R. Shared daughter: coupled to both parents. Run columns,
second half of the run: <E_par>, mean summed parent energy; swing, (max - min)/mean of it;
dissipation/driving, sum 2 gamma_d <E_d> over sum 2|gamma_p| <E_p>; n_active, daughters above 1e-3 of the
largest. stop: end (t = 1e13 s), runaway (summed energy 1e-2), wall (time limit).
"""


def lname(l: int, n: int) -> str:
    return f"({l},{n:+d})"


def rule(head: str, n_par: int, rc: bool) -> str:
    if "ranking" in head:
        r = "candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first"
        return r + ("; only pairs with an E_eq" if "only_eeq 1" in head else "")
    shared, direct = (int(re.search(rf"{k} (\d+)", head).group(1)) for k in ("shared", "direct"))
    d = (f"; plus, for each combination of the parents, its direct triplet of {'highest R_c' if rc else 'lowest E_eq'}"
         if direct else "")
    if shared:
        return "step k: the k daughters shared by both parents with the lowest larger E_th of their two pairs, each with both pairs" + d
    if n_par == 2:
        return "each parent's pairs ranked by E_th, the parents in turn" + d
    return "pairs ranked by E_th, lowest first; a pair that would bring N past the target skipped" + d


def selecting(lines: list[str]) -> list[list[str]]:
    rows = []
    for l in lines:
        if len(re.findall(MODE, l)) != 3 or "->" not in l:
            continue
        tok = re.sub(r"^#\s+(pair|harmonic|sum)?\s*", "", l).split()
        trip, nums = " ".join(tok[:5]), tok[5:]
        if tok[1] == "->":
            rows.append(["pair", trip, nums[0], nums[2], nums[3], nums[4], ""])
        else:
            kind = "harmonic" if tok[0] == tok[2] else "sum"
            old = len(nums) == 4                      # files before 2026-11-10 list R_c, not E_eq and mu
            rows.append([kind, trip, nums[0], nums[2], "", "" if old else nums[3],
                         f"(R_c {nums[3]})" if old else nums[-1]])
    return rows


def run_row(d: pathlib.Path) -> dict:
    modes = read_network(d / "network.data")[0]
    trip = [[int(x) for x in l.split()[:3]] for l in open(d / "network.data")
            if not l.startswith("#") and len(l.split()) == 5]
    par = [i for i, m in modes.items() if m["gen"] == 0]
    under: dict[int, set] = {}
    direct = []
    for t in trip:
        g = [modes[i]["gen"] for i in t]
        if g.count(0) == 2:
            x, y, c = sorted(t, key=lambda i: modes[i]["gen"])
            direct.append(f"{lname(modes[x]['l'], modes[x]['n'])} + {lname(modes[y]['l'], modes[y]['n'])}"
                          f" -> {lname(modes[c]['l'], modes[c]['n'])}")
        elif g.count(0) == 1:
            p = t[g.index(0)]
            for i in t:
                if i != p:
                    under.setdefault(i, set()).add(p)
    shared = [lname(modes[i]["l"], modes[i]["n"]) for i, ps in under.items() if len(ps) > 1]
    s = pd.read_csv(d / "summary.csv").iloc[0]
    e = d / "energies.csv"                      # not in git: absent on a fresh checkout
    t_stop = float(open(e).readlines()[-1].split(",")[0]) if e.exists() else math.nan
    return dict(run=d.name, N=len(modes) - len(par), triplets=len(trip), shared=" ".join(shared),
                direct="; ".join(direct), stop=s["stop"], t_stop=t_stop, E_par=s["E_parents"], swing=s["swing"],
                ratio=s["dissipation"] / s["driving"], n_active=s["n_active"])


def same(a: dict, b: dict) -> bool:
    return (a["stop"] == b["stop"] and a["n_active"] == b["n_active"] and abs(a["swing"] - b["swing"]) < 0.01
            and abs(a["E_par"] / b["E_par"] - 1) < 1e-3)


def ladder(lad: pathlib.Path) -> list[str]:
    runs = sorted((d for d in lad.iterdir() if (d / "summary.csv").exists()),
                  key=lambda d: int(re.match(r"d(\d+)", d.name).group(1)))
    if not runs:
        return []
    lines = open(runs[-1] / "network.data").read().splitlines()
    head = next(l for l in lines if l.startswith("# network_build"))
    modes = read_network(runs[-1] / "network.data")[0]
    par = [m for m in modes.values() if m["gen"] == 0]
    out = [f"### {lad.parent.name}/{lad.name}", ""]
    out += ["Parents: " + "; ".join(f"{lname(m['l'], m['n'])} f {m['w'] / (2 * math.pi) * 86400:.3f} c/d,"
                                    f" gamma {m['g']:+.3e} s^-1" for m in par), ""]
    out += [f"Built by: {rule(head, len(par), any('R_c' in l for l in lines))}.", ""]
    out += ["Selecting triplets of the largest network, in the order of its header:", "",
            "| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |",
            "|---|---|---|---|---|---|---|"]
    out += ["| " + " | ".join(r) + " |" for r in selecting(lines)]
    rows = [run_row(d) for d in runs]
    out += ["", "| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star]"
                " | swing | dissipation/driving | n_active |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    out += [f"| {r['run']} | {r['N']} | {r['triplets']} | {r['shared']} | {r['direct']} | {r['stop']} |"
            f" {r['t_stop']:.3g} | {r['E_par']:.4e} | {r['swing']:.3g} | {r['ratio']:.4f} | {r['n_active']} |"
            for r in rows]
    k = len(rows) - 1
    while k > 0 and same(rows[k - 1], rows[-1]):
        k -= 1
    out += ["", f"Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = {rows[k]['N']} on."
            if k < len(rows) - 1 else "Same run from some N on: not reached.", ""]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", type=pathlib.Path, required=True)
    a = ap.parse_args()
    out = [f"# {a.model.name}: ladders and runs", "", "Written by scripts/network_index.py.", "", NOTATION]
    for lp in ("ladder_1p", "ladder_2p"):
        if (a.model / lp).is_dir():
            out += [f"## {lp}", ""]
            for lad in sorted(p for p in (a.model / lp).iterdir() if p.is_dir()):
                out += ladder(lad)
    (a.model / "networks.md").write_text("\n".join(out) + "\n")
    print(f"-> {a.model}/networks.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
