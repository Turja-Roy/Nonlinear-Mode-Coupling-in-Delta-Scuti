"""Figures and tables for atol_bench.tex from out/atol_bench/{runs,plots}.

    python3 out/atol_bench/report/make_report.py
    cd out/atol_bench/report && latexmk -pdf atol_bench.tex
"""
import csv, pathlib, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle

R = pathlib.Path(__file__).resolve().parent
A, REPO = R.parent, R.parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from network_plot import read_network  # noqa: E402

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 6.5,
                     "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.02})

SOLVERS = {"rk8pd": "rk8pd", "bdf": "BDF, dense analytic Jacobian", "krylov_full": "BDF-Krylov, band full",
           "krylov_b3": "BDF-Krylov, band 3", "adams": "Adams, fixed point"}
ATOLS = ["1e-40", "1e-150", "1e-300", "0"]
RTOLS = ["1e-8", "1e-10", "1e-12"]
ATEX = {"1e-40": "10^{-40}", "1e-150": "10^{-150}", "1e-300": "10^{-300}", "0": "0",
        "1e-8": "10^{-8}", "1e-10": "10^{-10}", "1e-12": "10^{-12}"}
PH = "#efede8"
NOTES = {
    "1p2d_M1.7_eth": r"The second burst disagreed between rk8pd and BDF in the earlier sweep.",
    "1p2d_M2.0_dgamma": r"Deep dips between bursts.",
    "1p4d_M1.7_eth": r"Runaway near $3\times10^{11}$\,s.",
    "2p_M2.0T_eth": r"Daughter (14,$-48$) decays to $E\sim10^{-3600}$ in exact arithmetic; its regrowth "
                    r"time is set by each run's floor.",
    "2p_M2.0_eeq": r"Deep dips between bursts.",
}


def tex(s):
    return s.replace("_", r"\_")


def sci(t, p=2):
    m, e = f"{t:.{p}e}".split("e")
    return rf"${m}\times10^{{{int(e)}}}$"


def cnt(v):
    return "--" if not v else f"{v / 1e6:.2f}\\,M" if v >= 1e6 else f"{v / 1e3:.1f}\\,k" if v >= 1e3 else str(v)


def network(net):
    path = A / "nets" / f"{net}.data"
    head = open(path).readlines()[1].split()
    kv = dict(zip(head[4::2], head[5::2]))
    modes, pairs = [], []
    for line in open(path):
        if m := re.match(r"#\s+(\(\d+,[+-]\d+\)) -> (\(\d+,[+-]\d+\)) \+ (\(\d+,[+-]\d+\))\s+(.*)", line):
            pairs.append((*m.groups()[:3], *m[4].split()))
        f = line.split()
        if not line.startswith("#") and len(f) == 7:
            modes.append(dict(id=int(f[0]), gen=int(f[1]), lab=f"({int(float(f[3]))},{int(float(f[2])):+d})",
                              w=float(f[5]), g=float(f[6])))
    _, trip, _ = read_network(path)
    return dict(model=head[2].split("/")[-1], kv=kv, modes=modes, pairs=pairs, trip=trip)


def load(d):
    if not d.exists():
        return None
    if not (d / "summary.csv").exists():
        return "stalled"
    s = next(csv.DictReader(open(d / "summary.csv")))
    i = lambda k: int(float(s[k])) if s.get(k) not in (None, "") else None
    return dict(stop=s["stop"], tstop=float(s["t_end"]), sec=float(s["seconds"]), steps=i("steps"),
                rej=i("failed_steps"), rhs=i("rhs_calls"), jac=i("jac_evals"), nl=i("nonlin_iters"),
                lin=i("lin_iters"))


def fin(x):
    return isinstance(x, dict) and x["stop"] in ("end", "runaway")


def stats(x):
    if not isinstance(x, dict):
        return "--"
    s = rf"{cnt(x['steps'])} steps, {x['sec']:.1f}\,s"
    return s if x["stop"] == "end" else s + f", {x['stop']}"


def grid(rows, cols, png, runs):
    """rows x cols of network_plot.py images; a framed placeholder of the same size where none exists."""
    w, h = r"0.21\textwidth", r"0.1323\textwidth"            # network_plot.py's 1217 x 767 px
    out = [r"\begin{tabular}{@{}c@{\hspace{2pt}}" + r"c@{\hspace{2pt}}" * (len(cols) - 1) + r"c@{}}",
           " & " + " & ".join(rf"\small atol ${ATEX[c]}$" for c in cols) + r" \\[1pt]"]
    for rk, rl in rows:
        cells = [rf"\rotatebox{{90}}{{\parbox{{{h}}}{{\centering\scriptsize\bfseries {rl}}}}}"]
        for c in cols:
            p, x = png(rk, c), runs[(rk, c)]
            if p.exists():
                cells.append(rf"\includegraphics[width={w},height={h}]{{{p.relative_to(R, walk_up=True).as_posix()}}}")
            else:
                msg = "not run" if x is None else r"stalled\\[1pt]killed by the timeout,\\ no output written"
                cells.append(rf"{{\fboxsep0pt\fbox{{\parbox[b][\dimexpr{h}-2\fboxrule][c]{{\dimexpr{w}-2\fboxrule}}"
                             rf"{{\centering\scriptsize {msg}}}}}}}")
        out.append(" & ".join(cells) + r" \\[-1pt]")
        out.append(" & " + " & ".join(rf"\parbox[t]{{{w}}}{{\centering\scriptsize {stats(runs[(rk, c)])}}}" for c in cols) + r" \\[2pt]")
    out.append(r"\end{tabular}")
    return "\n".join(out)


def table(rows, first="integrator"):
    """rows: [(label, atol, run | 'stalled' | None)]."""
    out = [r"\begin{tabular}{@{}llllrrrrrrr@{}}", r"\toprule",
           rf"{first} & atol & stop & $t_{{\rm stop}}$ [s] & steps & rejected & RHS & Jac. & NL it. & lin. it. "
           r"& wall [s] \\", r"\midrule"]
    last = None
    for lab, al, x in rows:
        if last is not None and lab != last:
            out.append(r"\addlinespace[2pt]")
        name = lab if lab != last else ""
        last = lab
        if not isinstance(x, dict):
            out.append(f"{name} & ${ATEX[al]}$ & {'stalled' if x else 'not run'} & \\multicolumn{{8}}{{l}}"
                       f"{{{'no output written' if x else ''}}} \\\\")
            continue
        out.append(f"{name} & ${ATEX[al]}$ & {x['stop']} & {sci(x['tstop'])} & {cnt(x['steps'])} & "
                   f"{cnt(x['rej'])} & {cnt(x['rhs'])} & {cnt(x['jac'])} & {cnt(x['nl'])} & {cnt(x['lin'])} & "
                   f"{x['sec']:.1f} \\\\")
    out += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(out)


def header(net, nw):
    kv = nw["kv"]
    mn = lambda x: x.replace("-", "$-$")
    modes = "\n".join(f"{m['id']} & {m['gen']} & {mn(m['lab'])} & {sci(m['w'], 3)} & {sci(m['g'], 3)} \\\\"
                      for m in nw["modes"])
    pairs = "\\\\\n".join(f"{mn(a)} $\\to$ {mn(b)} + {mn(c)}: $\\Delta/\\omega_a={sci(float(r))[1:-1]}$, "
                          f"$\\kappa={float(k):+.3f}$, $E_{{\\rm th}}={sci(float(et))[1:-1]}$, "
                          + ("no $E_{\\rm eq}$" if ee == "none" else f"$E_{{\\rm eq}}={sci(float(ee))[1:-1]}$")
                          for a, b, c, _, r, k, et, ee in nw["pairs"])
    return rf"""\clearpage
\subsection{{\texttt{{{tex(net)}}}}}\label{{net:{net}}}
\noindent\begin{{minipage}}[t]{{0.43\textwidth}}\vspace{{0pt}}\footnotesize
\begin{{tabular}}{{@{{}}rrlrr@{{}}}}
\toprule
id & gen & $(\ell,n)$ & $\omega$ [rad\,s$^{{-1}}$] & $\gamma$ [s$^{{-1}}$] \\
\midrule
{modes}
\bottomrule
\end{{tabular}}
\end{{minipage}}\hfill
\begin{{minipage}}[t]{{0.55\textwidth}}\vspace{{0pt}}\footnotesize
\textbf{{Model}} \texttt{{{tex(nw['model'])}}}, rule {kv.get('rule')}, $N={kv.get('N')}$, cut {kv.get('cut')},
closure {float(kv.get('closure', 0)):g}; {len(nw['modes'])} modes, {len(nw['trip'])} triplet(s).\\[2pt]
\textbf{{Selecting pairs}}\\ {pairs}\\[2pt]
{NOTES.get(net, '')}
\end{{minipage}}
"""


def maps(path, nets, cells):
    """Steps and wall time for every network x (integrator, atol) at rtol 1e-8."""
    cols = [(s, a) for s in SOLVERS for a in ATOLS]
    fig, axs = plt.subplots(2, 1, figsize=(7.1, 5.2), layout="constrained")
    specs = [("accepted steps", [0, 1e6, 3e6, 1e7, 3e7, 1e8, 1e12], "Purples", "steps",
              lambda v: f"{v / 1e6:.1f}M" if v < 1e7 else f"{v / 1e6:.0f}M",
              ["<1M", "1–3M", "3–10M", "10–30M", "30–100M", ">100M"]),
             ("wall time [s]", [0, 3, 10, 30, 100, 300, 1e9], "Blues", "sec",
              lambda v: f"{v:.0f}" if v >= 10 else f"{v:.1f}", ["<3", "3–10", "10–30", "30–100", "100–300", ">300"])]
    for ax, (name, bounds, cm, key, fmt, ticks) in zip(axs, specs):
        k = len(bounds) - 1
        cmap = ListedColormap(plt.get_cmap(cm)(np.linspace(0.12, 0.85, k)))
        norm = BoundaryNorm(bounds, k)
        M = np.full((len(nets), len(cols)), np.nan)
        for i, net in enumerate(nets):
            for j, (s, a) in enumerate(cols):
                x = cells[net][(s, a)]
                if fin(x):
                    M[i, j] = x[key]
                    ax.text(j, i, fmt(x[key]), ha="center", va="center", fontsize=4.6,
                            color="white" if norm(x[key]) >= k - 2 else "#222")
                else:
                    ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fc=PH, ec="#c9c6bf", lw=0, hatch="////"))
                    ax.text(j, i, "stalled" if x == "stalled" else x["stop"], ha="center", va="center",
                            fontsize=4.3, color="#222")
        im = ax.imshow(M, cmap=cmap, norm=norm, aspect="auto")
        ax.set_xticks(range(len(cols)), [a for _, a in cols], fontsize=5, rotation=90)
        ax.set_yticks(range(len(nets)), nets, fontsize=5.5)
        for j in range(4, len(cols), 4):
            ax.axvline(j - 0.5, color="white", lw=2.5)
        sec = ax.secondary_xaxis("top")
        sec.set_xticks([1.5 + 4 * g for g in range(len(SOLVERS))], list(SOLVERS.values()), fontsize=6)
        sec.tick_params(length=0)
        ax.tick_params(length=0)
        cb = fig.colorbar(im, ax=ax, ticks=[(bounds[q] + bounds[q + 1]) / 2 for q in range(k)],
                          fraction=0.025, pad=0.01)
        cb.ax.set_yticklabels(ticks, fontsize=5.5)
        cb.ax.tick_params(length=0)
        cb.set_label(name, fontsize=6.5)
    fig.savefig(path)
    plt.close(fig)


def main():
    (R / "figs").mkdir(exist_ok=True)
    (R / "gen").mkdir(exist_ok=True)
    nets = sorted((d.name for d in (A / "runs/rk8pd").iterdir()), key=lambda n: (n.split("_")[0][:2], n))
    cells, pages, appendix = {}, [], []
    for net in nets:
        nw = network(net)
        runs = {(s, a): load(A / "runs" / s / net / f"r1e-8_a{a}") for s in SOLVERS for a in ATOLS}
        cells[net] = runs
        g = grid(list(SOLVERS.items()), ATOLS, lambda s, a: A / "plots" / s / net / f"r1e-8_a{a}.png", runs)
        t = table([(SOLVERS[s], a, runs[(s, a)]) for s in SOLVERS for a in ATOLS])
        pages.append(header(net, nw) + rf"""
\par\smallskip{{\centering{g}\par}}
\par\noindent{{\scriptsize\setlength\tabcolsep{{3.5pt}}\renewcommand\arraystretch{{0.92}}
{t}}}
""")
        runs = {(r, a): load(A / "runs/rk8pd" / net / f"r{r}_a{a}") for r in RTOLS for a in ATOLS}
        g = grid([(r, f"rtol ${ATEX[r]}$") for r in RTOLS], ATOLS,
                 lambda r, a: A / "plots/rk8pd" / net / f"r{r}_a{a}.png", runs)
        t = table([(f"${ATEX[r]}$", a, runs[(r, a)]) for r in RTOLS for a in ATOLS], first="rtol")
        appendix.append(rf"""\clearpage
\subsection{{\texttt{{{tex(net)}}}: rk8pd, rtol $\times$ atol}}
\noindent See \S\ref{{net:{net}}} for the network. {NOTES.get(net, '')}
\par\smallskip{{\centering{g}\par}}
\par\noindent{{\scriptsize\setlength\tabcolsep{{3.5pt}}\renewcommand\arraystretch{{0.92}}
{t}}}
""")

    maps(R / "figs/maps.pdf", nets, cells)
    (R / "gen/networks.tex").write_text("\n".join(pages))
    (R / "gen/appendix.tex").write_text("\n".join(appendix))

    lines = [r"\begin{tabular}{@{}llrrrl@{}}", r"\toprule",
             r"network & model & modes & parents & triplets & rule \\", r"\midrule"]
    for net in nets:
        nw = network(net)
        lines.append(f"\\hyperref[net:{net}]{{\\texttt{{{tex(net)}}}}} & \\texttt{{{tex(nw['model'])}}} & "
                     f"{len(nw['modes'])} & {sum(m['gen'] == 0 for m in nw['modes'])} & {len(nw['trip'])} & "
                     f"{nw['kv'].get('rule')} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (R / "gen/netlist.tex").write_text("\n".join(lines))

    lines = [r"\begin{tabular}{@{}llrrrrr@{}}", r"\toprule",
             r"integrator & atol & stopped early & steps / rk8pd & wall / rk8pd & $\sum$ steps & $\sum$ wall [s] \\",
             rf" & & (of {len(nets)}) & median & median & & \\", r"\midrule"]
    for s in SOLVERS:
        for a in ATOLS:
            xs = [cells[n][(s, a)] for n in nets]
            both = [(x, cells[n][("rk8pd", a)]) for n, x in zip(nets, xs) if fin(x) and fin(cells[n][("rk8pd", a)])]
            done = [x for x in xs if fin(x)]
            tot = "" if len(done) == len(nets) else r"$^*$"
            lines.append(f"{SOLVERS[s] if a == ATOLS[0] else ''} & ${ATEX[a]}$ & {len(nets) - len(done)} & "
                         f"{np.median([x['steps'] / x0['steps'] for x, x0 in both]):.2f} & "
                         f"{np.median([x['sec'] / x0['sec'] for x, x0 in both]):.2f} & "
                         f"{cnt(sum(x['steps'] for x in done))}{tot} & {sum(x['sec'] for x in done):.0f}{tot} \\\\")
        lines.append(r"\addlinespace[2pt]")
    lines[-1:] = [r"\bottomrule", r"\end{tabular}"]
    (R / "gen/summary.tex").write_text("\n".join(lines))
    print((R / "gen/summary.tex").read_text())


if __name__ == "__main__":
    main()
