#!/usr/bin/env python3
"""Figures for the progress report (Reading/Papers/Study/figs/).

Every panel comes from data already on disk; the only computation is the
cumulative kappa of two representative triplets (milliseconds each).
"""

from __future__ import annotations

import functools
import pathlib
import sys
from typing import NamedTuple

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from coupling.kappa import kappa_abc
from coupling.modes import gamma_turb, load_model
from coupling.observables import classify_frame, daughter_index

CD = 7.27220522e-5  # rad/s per cycle/day
DETUNING_CUT = 0.15  # in units of sqrt(GM/R^3), as in three_mode_search.py
MU_MIN = None  # lower y limit of every mu panel; None keeps the full range

# argparse defaults; main() sets the rest.
MODEL = pathlib.Path("models/dsct_M2.0")
MODEL_ROOT = pathlib.Path("models")  # for cross-model figures; main() resets it
OUT = pathlib.Path("out")
TAG = FIGS = bg = efs = None
SOURCE = "narrow"   # which coupling table _channel_table reads
L_MAX_CUT = None    # keep rows with max(l_b, l_c) <= this

# Okabe-Ito, fixed assignment order (colorblind-safe)
BLUE, ORANGE, GREEN, VERM, PURPLE, SKY = (
    "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9")

plt.rcParams.update({
    "figure.dpi": 150, "savefig.bbox": "tight",
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 11,
    "lines.linewidth": 1.6,
})


def save(fig, name):
    # A ladder rung has to say which rung it is: these figures are meant to be
    # opened side by side, where the directory name is no longer visible. The
    # restriction matters too -- the sweep ran --parents-only, so slot a is
    # always driven and direct-sum, inactive and all-damped cannot appear.
    if SOURCE == "lconv":
        cut = (r"full net, $l \leq 25$" if L_MAX_CUT is None
               else rf"$l_{{\rm max}} = {L_MAX_CUT}$")
        fig.text(0.995, -0.015,
                 f"wide sweep, {cut};  driven sum mode only,  "
                 r"$m = (0,0,0)$,  $\gamma_{\rm rad}+\gamma_{\rm turb}$",
                 ha="right", va="top", fontsize=7, color="0.45")
    for ext in ("pdf", "png"):
        fig.savefig(FIGS / f"{name}.{ext}")
    plt.close(fig)
    print(f"  {name}")


def load_h5(path):
    out = {}
    with h5py.File(path, "r") as f:
        out.update(f.attrs)
        for k in f:
            v = f[k][...]
            out[k] = v["re"] + 1j * v["im"] if v.dtype.names == ("re", "im") else v
    return out




def _core_boundary(x, N2):
    """Outer edge of the convective core: first x where N^2 turns positive."""
    pos = np.flatnonzero(N2 > 0)
    return float(x[pos[0]]) if pos.size else float("nan")


def _daughter_band(drv_lo, drv_hi):
    """Parametric-daughter band in c/d.

    Prefer the gnet run's own frequency span -- pass 2 is pure gnet, so its
    summary *is* the band.  Without it, fall back to half the driven band,
    which is the omega_c/2 resonance condition with the direct daughter taken
    at the parent frequencies.
    """
    path = MODEL / "gyre" / "summary_nad_wide_hi.h5"
    if path.exists():
        f = np.real(load_h5(path)["freq"])
        return float(f.min()), float(f.max())
    return drv_lo / 2.0, drv_hi / 2.0


def _brunt2():
    """(x, r, N^2) outside the innermost points, where g -> 0 makes N^2 noisy."""
    m = bg.x > 1e-3
    x, r = bg.x[m], bg.r[m]
    dlnP_dr = -bg.rho[m] * bg.g[m] / bg.P[m]
    dlnrho_dr = bg.dlnrho_dlnr[m] / r
    return x, r, -bg.g[m] * (dlnrho_dr - dlnP_dr / bg.Gamma_1[m])


def _x_core():
    x, _, N2 = _brunt2()
    return _core_boundary(x, N2)


# ---------------------------------------------------------------- propagation
def fig_propagation():
    x, r, N2 = _brunt2()
    m = bg.x > 1e-3
    cs2 = bg.Gamma_1[m] * bg.P[m] / bg.rho[m]
    to_cd = 86400 / (2 * np.pi)

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    N = np.sqrt(np.clip(N2, 0, None)) * to_cd
    ax.plot(x, np.where(N > 0, N, np.nan), color=BLUE, label=r"$N$")
    for l, ls in ((1, ":"), (2, "--"), (3, "-")):
        Sl = np.sqrt(l * (l + 1) * cs2) / r * to_cd
        ax.plot(x, Sl, color=ORANGE, ls=ls, lw=1.2, label=rf"$S_{{{l}}}$")
    # Bands and core boundary come from this model, not from dsct_M2.0's.
    drv_lo, drv_hi = (v / CD for v in _gamma_sets()[3])
    dtr_lo, dtr_hi = _daughter_band(drv_lo, drv_hi)
    x_core = _core_boundary(x, N2)

    ax.axhspan(drv_lo, drv_hi, color=GREEN, alpha=0.12, lw=0)
    ax.text(0.55, np.sqrt(drv_lo * drv_hi), "linearly driven band",
            color=GREEN, fontsize=9)
    ax.axhspan(dtr_lo, dtr_hi, color=VERM, alpha=0.12, lw=0)
    ax.text(0.55, np.sqrt(dtr_lo * dtr_hi) * 0.85,
            r"parametric daughter band ($\omega_c/2$)", color=VERM, fontsize=9)
    ax.axvline(x_core, color="0.4", ls=":", lw=1)
    ax.text(x_core + 0.007, 0.28, "conv. core\nboundary", fontsize=8, color="0.35")
    ax.set(xlabel=r"$x = r/R$", ylabel="frequency (c/d)", yscale="log",
           xlim=(0, 1.0), ylim=(0.2, 300),
           title=rf"Propagation diagram, $M = {bg.M / 1.988409870698051e33:.2f}"
                 rf"\,M_\odot$ $\delta$ Sct model")
    ax.legend(loc="upper right", ncols=2, framealpha=0.9)
    save(fig, "fig_propagation")


# ------------------------------------------------------------------- damping
def _gamma_sets():
    """Point sets for the damping/driving figures, in MW23 Fig. 4 units.

    GYRE reports the dimensionless eigenvalue omega = sigma / sqrt(GM/R^3),
    so gamma[s^-1] = -Im(omega) * sqrt(GM/R^3), and freq (c/d) -> rad/s is
    the usual 2*pi/86400.  Unlike MW23 Fig. 4 we keep the linearly unstable
    modes on the plot; they drop their radiative point but keep their
    turbulent one, which is positive for every mode.

    The turbulent set comes from the narrow-net eigenfunctions, so it is only
    the l <= 3 modes with a detail dump -- the summary files carry no
    displacement.

    Returns (hi_l, lo_damped, lo_driven, band, turb, comb): (omega, |gamma|)
    pairs, plus `comb` = (damped, driven, n_flip) for the low-l total
    gamma_rad + gamma_turb. `band` stays radiative, so the shading matches
    across figures.
    """
    G = 6.67430e-8  # cgs
    nar = load_h5(MODEL / "gyre/summary_nad.h5")
    # The wide passes are optional: not every model has had them run.
    wide = [load_h5(p) for p in (MODEL / "gyre" / n for n in
                                 ("summary_nad_wide.h5", "summary_nad_wide_hi.h5"))
            if p.exists()]
    w_dyn = np.sqrt(G * nar["M_star"][0] / nar["R_star"][0] ** 3)  # s^-1

    def unpack(d):
        return (np.real(d["freq"]) * CD, -np.imag(d["omega"]) * w_dyn,
                d["l"], d["n_pg"])

    cols = [unpack(nar)] + [unpack(d) for d in wide]
    w, g, l, npg = (np.concatenate([c[i] for c in cols]) for i in range(4))

    hi = (l > 3) & (g > 0)
    hi_l = (w[hi], g[hi])

    # low-l set = narrow net plus the deeper l <= 3 g-modes of the wide net,
    # deduplicated on (l, n_pg), so no damped low-l mode is hidden by the
    # l > 3 cut above.
    lo = l <= 3
    _, uniq = np.unique(np.stack([l[lo], npg[lo]]), axis=1, return_index=True)
    w, g, l_lo, npg_lo = w[lo][uniq], g[lo][uniq], l[lo][uniq], npg[lo][uniq]
    dam, drv = g > 0, g < 0
    band = (w[drv].min(), w[drv].max())

    # One pass: gamma_turb is a spline plus a quadrature per mode, and both the
    # turbulent set and the combined one need it.
    gt = {k: gamma_turb(ef) for k, ef in efs.items()}
    turb = np.array([[efs[k].omega, v] for k, v in gt.items()]).T
    # Wide-pass modes have no detail dump and so no turbulent rate; they keep
    # their radiative value rather than dropping off the plot.
    tur = np.array([gt.get((int(li), int(ni)), 0.0)
                    for li, ni in zip(l_lo, npg_lo)])
    tot = g + tur
    n_flip = int(((g < 0) & (tot > 0)).sum())
    return (hi_l, (w[dam], g[dam]), (w[drv], -g[drv]), band, (turb[0], turb[1]),
            ((w[tot > 0], tot[tot > 0]), (w[tot < 0], -tot[tot < 0]), n_flip))


LBL_HI = r"damped, $4 \leq l \leq 25$ (daughter net)"
LBL_DAM = r"damped, $l \leq 3$"
LBL_DRV = r"driven ($\gamma < 0$), $l \leq 3$"
LBL_TRB = r"turbulent, $l \leq 3$"
LBL_TOT_DAM = r"damped, $l \leq 3$  ($\gamma_{\rm rad} + \gamma_{\rm turb}$)"
LBL_TOT_DRV = r"driven, $l \leq 3$  ($\gamma_{\rm tot} < 0$)"


def _gamma_axes(ax, band, top_axis=True):
    ax.axvspan(*band, color=VERM, alpha=0.07, lw=0)
    ax.set(xscale="log", yscale="log")
    if top_axis:
        sec = ax.secondary_xaxis("top", functions=(lambda x: x / CD,
                                                   lambda x: x * CD))
        sec.set_xlabel("frequency (c/d)", fontsize=9)


def fig_gamma():
    """Low-l total gamma_rad + gamma_turb; fig_gamma_panels keeps the split.
    The l > 3 net has no detail dump, so its points stay radiative."""
    hi_l, _dam, _drv, band, _trb, (tot_dam, tot_drv, n_flip) = _gamma_sets()
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.scatter(*hi_l, s=4, color="0.75", lw=0, rasterized=True,
               label=LBL_HI + ", radiative only")
    ax.scatter(*tot_dam, s=14, color=BLUE, lw=0, label=LBL_TOT_DAM)
    ax.scatter(*tot_drv, s=22, color=VERM, marker="D", lw=0, label=LBL_TOT_DRV)
    _gamma_axes(ax, band)
    ax.set(xlabel=r"mode frequency $\omega$  $[\mathrm{rad\,s^{-1}}]$",
           ylabel=r"total rate $|\gamma_{\rm rad} + \gamma_{\rm turb}|$"
                  r"  $[\mathrm{s^{-1}}]$",
           title=r"Total damping rate, radiative plus turbulent")
    if n_flip:
        ax.text(0.02, 0.97, f"{n_flip} radiatively driven mode"
                            f"{'s' if n_flip > 1 else ''} net damped once\n"
                            r"$\gamma_{\rm turb}$ is included",
                transform=ax.transAxes, va="top", fontsize=8, color=VERM)
    ax.legend(loc="lower right", framealpha=0.9, fontsize=8)
    save(fig, "fig_gamma")


def fig_gamma_panels():
    """Same data as fig_gamma, split by population.

    The point of panel (b) is that the *radiative* points are empty inside the
    shaded driven band: there is no radiatively damped l <= 3 mode there for
    the diamonds to hide. The turbulent points do fill the band, since
    turbulent viscosity damps every mode including the driven ones; they sit
    well below the radiative points, which is MW23's gamma_tot ~ gamma_rad.
    """
    hi_l, dam, drv, band, trb, _comb = _gamma_sets()
    n_in = lambda s: int(((s[0] >= band[0]) & (s[0] <= band[1])).sum())

    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.0), sharex=True, sharey=True)
    panels = [
        (axes[0, 0], "(a) daughter net, $4 \\leq l \\leq 25$ (all damped)",
         [(hi_l, dict(s=5, color="0.6", lw=0, rasterized=True), LBL_HI)],
         f"{n_in(hi_l)} in band", "lower right"),
        (axes[0, 1], "(b) damped, $l \\leq 3$: radiative vs turbulent",
         [(dam, dict(s=18, color=BLUE, lw=0), LBL_DAM),
          (trb, dict(s=12, color=PURPLE, marker="^", lw=0), LBL_TRB)],
         f"{n_in(dam)} radiative in band", "upper left"),
        (axes[1, 0], "(c) driven, $l \\leq 3$",
         [(drv, dict(s=26, color=VERM, marker="D", lw=0), LBL_DRV)],
         f"{n_in(drv)} in band", "upper left"),
        (axes[1, 1], "(d) all together",
         [(hi_l, dict(s=4, color="0.75", lw=0, rasterized=True), LBL_HI),
          (dam, dict(s=14, color=BLUE, lw=0), LBL_DAM),
          (drv, dict(s=22, color=VERM, marker="D", lw=0), LBL_DRV),
          (trb, dict(s=10, color=PURPLE, marker="^", lw=0), LBL_TRB)],
         None, "upper left"),
    ]
    for ax, title, layers, note, legloc in panels:
        for pts, kw, lbl in layers:
            ax.scatter(*pts, label=lbl, **kw)
        _gamma_axes(ax, band, top_axis=ax in (axes[0, 0], axes[0, 1]))
        ax.set_title(title, fontsize=10)
        ax.legend(loc=legloc, framealpha=0.9, fontsize=7)
        if note:
            ax.text(np.sqrt(band[0] * band[1]), 3e-5, note, ha="center",
                    fontsize=9, color=VERM if note.startswith("0") else "0.3",
                    fontweight="bold" if note.startswith("0") else "normal")
    for ax in axes[1]:
        ax.set_xlabel(r"mode frequency $\omega$  $[\mathrm{rad\,s^{-1}}]$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$|\gamma|$  $[\mathrm{s^{-1}}]$")
    fig.suptitle("Damping and driving rates by population "
                 "(shaded = driven band, %.2f--%.2f c/d)"
                 % (band[0] / CD, band[1] / CD), y=1.0, fontsize=11)
    save(fig, "fig_gamma_panels")


# -------------------------------------------------------------- kappa_cum
@functools.cache
def _kappa_table():
    """kappa_<tag>.csv at m = 0; kappa depends on m only through the angular
    scalars, so one m combination stands for all."""
    kap = pd.read_csv(OUT / f"kappa_{TAG}.csv")
    return kap[(kap.m_a == 0) & (kap.m_b == 0) & (kap.m_c == 0)]


def _best_rows(by: str = "mu"):
    """The all-p and all-g triplet maximising `by`, in that order. Mixed-sign
    triplets are dropped; they blur the surface-versus-core contrast.

    "mu" and "abs_kappa" select different triplets, since mu also carries the
    detuning and the damping.
    """
    kap = _kappa_table()
    allpos = (kap.n_a > 0) & (kap.n_b > 0) & (kap.n_c > 0)
    allneg = (kap.n_a < 0) & (kap.n_b < 0) & (kap.n_c < 0)
    return [kap[m].loc[kap[m][by].idxmax()] for m in (allpos, allneg)]


def _triplet_keys(r):
    return [(int(r.l_a), int(r.n_a)), (int(r.l_b), int(r.n_b)),
            (int(r.l_c), int(r.n_c))]


def _sci(v, digits=2):
    """1.4e+04 -> 1.4 \\times 10^{4}, for math mode. Moderate values pass
    through unchanged."""
    if 1e-2 <= abs(v) < 1e4:
        return f"{v:.{digits + 1}g}"
    mant, exp = f"{v:.{digits - 1}e}".split("e")
    return rf"{mant} \times 10^{{{int(exp)}}}"


def fig_kappa_cum():
    rows = _best_rows()
    x_core = _x_core()

    fig, axes = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True)
    for ax, r, color, tag in zip(axes, rows, (BLUE, VERM),
                                 ("p-mode triplet", "g-mode triplet")):
        keys = _triplet_keys(r)
        res = kappa_abc(*(efs[k] for k in keys), (0, 0, 0))
        ax.plot(res.r / bg.R, res.cumulative, color=color)
        ax.axhline(0, color="0.6", lw=0.7)
        ax.axvline(x_core, color="0.4", ls=":", lw=1)
        lbl = ", ".join(f"$({l},{n:+d})$" for l, n in keys)
        ax.set_title(f"{tag}: {lbl},  $\\kappa = {res.kappa:.3g}$", fontsize=10)
        ax.set_ylabel(r"$\kappa(<r)$")
    axes[0].text(x_core, axes[0].get_ylim()[1] * 0.75, " conv. core boundary",
                 fontsize=8, color="0.35")
    axes[1].set_xlabel(r"$x = r/R$")
    fig.suptitle(r"Cumulative coupling integral: $p$ builds at the surface, "
                 r"$g$ at the core boundary", y=0.99, fontsize=11)
    save(fig, "fig_kappa_cum")


# ------------------------------------------------ xi_r + kappa(<r), MW23 2 & 3
CRITERIA = (("mu", r"largest $\mu$"), ("abs_kappa", r"largest $|\kappa_{abc}|$"))


def _xi_kappa_column(top, bot, row, depth: bool, color: str):
    """One (|xi_r|, kappa(<r)) column for one triplet."""
    keys = _triplet_keys(row)
    res = kappa_abc(*(efs[k] for k in keys), (0, 0, 0))
    l, n = max(keys, key=lambda k: abs(k[1]))
    ef = efs[(l, n)]
    R, x_core = bg.R, _x_core()
    abscissa = (lambda r: R - r) if depth else (lambda r: r)

    cut = ef.bg.x <= 1.0
    for ax, u, v in ((top, abscissa(ef.bg.r[cut]), np.abs(ef.xi_r[cut])),
                     (bot, abscissa(res.r), res.cumulative)):
        keep = u > 0  # the origin of a log abscissa: r = R for depth, r = 0 for radius
        ax.plot(u[keep], v[keep], color=color)
        ax.axvline(abscissa(x_core * R), color="0.4", ls=":", lw=1)
    bot.axhline(0, color="0.6", lw=0.7)
    bot.axhline(res.kappa, color="0.6", ls="--", lw=0.7)

    # Label goes opposite the eigenfunction's peak.
    y_lbl, va = (0.97, "top") if depth else (0.03, "bottom")
    top.text(0.02, y_lbl, rf"$l = {l}$, $n_{{pg}} = {n:+d}$ "
             + ("$p$" if n > 0 else "$g$") + "-mode",
             transform=top.transAxes, va=va, fontsize=9, color=color)
    return res


def _xi_kappa(depth: bool, name: str):
    """MW23 Fig. 2 (depth=True) and Fig. 3 (depth=False), one column per entry
    of CRITERIA.

    Top row is |xi_r| of the triplet's largest-|n_pg| member -- a log axis
    cannot show the nodes' sign flips -- bottom row kappa_abc(<r) of the whole
    triplet. `depth` puts the surface on the right, as MW23 do for p-modes.
    kappa(<r) changes sign, so only the top row is log-log, and the two kappa
    panels keep independent ordinates.
    """
    R, x_core = bg.R, _x_core()
    abscissa = (lambda r: R - r) if depth else (lambda r: r)
    lo = 1e7 if depth else 0.5 * x_core * R
    xlim = (R, lo) if depth else (lo, R)
    colors = (BLUE, ORANGE) if depth else (VERM, PURPLE)
    rows = [_best_rows(by)[0 if depth else 1] for by, _ in CRITERIA]

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 6.8), sharex=True,
                             gridspec_kw=dict(hspace=0.08, wspace=0.2))
    axes[0, 1].sharey(axes[0, 0])
    axes[0, 1].tick_params(labelleft=False)

    for j, (row, color, (_, crit)) in enumerate(zip(rows, colors, CRITERIA)):
        top, bot = axes[0, j], axes[1, j]
        res = _xi_kappa_column(top, bot, row, depth, color)
        top.set(xscale="log", yscale="log", xlim=xlim)
        lbl = ", ".join(f"$({la},{na:+d})$" for la, na in _triplet_keys(row))
        top.set_title(f"{crit}:  {lbl}\n"
                      rf"$\kappa_{{abc}} = {res.kappa:.3g}$,"
                      rf"  $\mu = {_sci(row.mu)}$", fontsize=10)
        bot.set_xlabel(r"$R - r$  $[\mathrm{cm}]$" if depth
                       else r"$r$  $[\mathrm{cm}]$")

    axes[0, 0].set_ylabel(r"$|\xi_r(r)|$  $[\mathrm{cm}]$"
                          "\n" r"(normalised to $E = E_\star$)")
    for j in (0, 1):
        axes[1, j].set_ylabel(r"$\kappa_{abc}(<r)$")
    # Offset in points, not data units: the boundary sits near an axis edge.
    axes[0, 0].annotate("conv. core boundary", (abscissa(x_core * R), 0.5),
                        xycoords=("data", "axes fraction"),
                        textcoords="offset points", xytext=(5, 0), rotation=90,
                        ha="left", va="center", fontsize=8, color="0.35")
    save(fig, name)


def fig_xi_kappa_p():
    _xi_kappa(depth=True, name="fig_xi_kappa_p")


def fig_xi_kappa_g():
    _xi_kappa(depth=False, name="fig_xi_kappa_g")


# ----------------------------------------------------- combination amplitudes
def fig_mu():
    obs = _read_observables(["omega_a", "omega_b", "omega_c",
                             "mu_a", "mu_b", "mu_c", "mu_max"])
    mus = obs[["mu_a", "mu_b", "mu_c"]].to_numpy()
    ws = np.abs(obs[["omega_a", "omega_b", "omega_c"]].to_numpy()) / CD
    f_driven = ws[np.arange(len(obs)), mus.argmax(axis=1)]
    A = obs.mu_max.to_numpy() * 1e-12  # A_x = mu A_a A_b at A = 1e-6

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.scatter(f_driven, A, s=3, color=BLUE, alpha=0.25, lw=0, rasterized=True)
    ax.axhspan(2e-12, 8e-12, color=ORANGE, alpha=0.25, lw=0)
    ax.text(30, 1.2e-13, "ordinary combination\nfrequencies ($\\mu \\approx 4$)",
            fontsize=9, color=VERM)
    ax.annotate("", xy=(52, 4e-12), xytext=(48, 3e-13),
                arrowprops=dict(arrowstyle="->", color=VERM, lw=1))
    ax.set(yscale="log", xlabel="driven-mode frequency (c/d)", xlim=(0, 85),
           ylabel=r"$A_x = \mu\, A_a A_b$   (at $A_{a,b} = 10^{-6}$)",
           title=r"Predicted three-mode amplitudes, all 71142 triplets")
    save(fig, "fig_mu")


# ------------------------------------------------------------ E_th histogram
def fig_eth():
    obs = _read_observables(["gamma_a", "E_th_over_E_star"])
    eth = obs.E_th_over_E_star[(obs.gamma_a < 0) & (obs.E_th_over_E_star > 0)]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.hist(np.log10(eth), bins=60, color=BLUE, alpha=0.85)
    # ax.axvline(np.log10(eth.min()), color=VERM, ls="--", lw=1.4)
    # ax.text(np.log10(eth.min()) + 0.15, ax.get_ylim()[1] * 0.85,
    #         f"min $= {eth.min():.2g}$", color=VERM, fontsize=9)
    # ax.axvline(-12, color=GREEN, ls="--", lw=1.4)
    # ax.text(-11.7, ax.get_ylim()[1] * 0.45,
    #         "observed parent\nenergies $\\sim 10^{-12}$", color=GREEN, fontsize=9)
    # ax.set(xlabel=r"$\log_{10} (E_{\rm th}/E_\star)$", ylabel="triplets",
    #        title="Parametric thresholds, driven parents: floor sits "
    #              r"$1660\times$ above observed")
    ax.set_ylabel("triplets")
    ax.set_xlabel(r"$\log_{10} (E_{\rm th}/E_\star)$")
    save(fig, "fig_eth")


def fig_eth_pg():
    """fig_eth split by p/g make-up. Every row here is parametric: E_th is zero
    unless the sum mode is the lone parent (`observables.threshold_energy`), so
    the make-up is parent + its two daughters. The title reports the daughter
    census rather than assuming it -- where every daughter is a g-mode the
    make-up reads directly as the parent's type."""
    df = _channel_table()
    eth = df[(df.gamma_a < 0) & (df.E_th_over_E_star > 0)]
    comp = _pg_comp(eth[["n_a", "n_b", "n_c"]].to_numpy())
    x = np.log10(eth.E_th_over_E_star.to_numpy())
    dau = _pg_letters(eth[["n_b", "n_c"]].to_numpy()).ravel()
    n_g = int((dau == "g").sum())

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    _pg_stack(ax, x, comp, np.linspace(x.min(), x.max(), 60))
    ax.set(xlabel=r"$\log_{10} (E_{\rm th}/E_\star)$", ylabel="triplets",
           title=f"Parametric thresholds by $p$/$g$ make-up, {len(eth)} rows\n"
                 f"{n_g} of {len(dau)} daughters are $g$-modes")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    save(fig, "fig_eth_pg")


# ------------------------------------------------------ coupling channels
# (colour, alpha, legend label). fig_detuning reads the colour only.
CHANNEL_STYLE = {
    "direct-diff": (VERM, 0.5, "difference: $\\omega_a - \\omega_x$ drives a damped daughter"),
    "direct-sum": (BLUE, 0.5, "sum: two driven daughters drive the parent"),
    "parametric": (GREEN, 0.5, "parametric: driven parent decays into two daughters"),
}


# MW23 Fig. 5's triplet types. p-mode is n_pg > 0, g-mode n_pg < 0 -- their own
# convention (sec 3.1: "-20 <= n <= 20, corresponding to g- and p-modes").
# ponytail: sign(n_pg) only. ~20 modes in the observed band are genuinely mixed
# (n_p > 0 and n_g > 0), where the Eckart label is unreliable; n_p, n_g sit in
# summary_nad.h5 and can be read here with load_h5 if that ever matters.
def _pg_type(ns, idx, j, rows):
    """The five types, empty where they do not apply: no direct channel to name
    a single daughter, or an f-mode (n_pg = 0) in the triplet."""
    p = ns > 0
    n_par_p = p.sum(axis=1) - p[rows, j]  # p-modes among the two parents
    dau_p = p[rows, j]
    out = np.full(len(rows), "mixed", dtype=object)
    out[(n_par_p == 2) & dau_p] = "ppp"
    out[(n_par_p == 2) & ~dau_p] = "ppg"
    out[(n_par_p == 0) & ~dau_p] = "ggg"
    out[(n_par_p == 0) & dau_p] = "ggp"
    out[(ns == 0).any(axis=1)] = "f"  # f-mode: no p/g label
    out[idx < 0] = ""  # no direct channel at all; this wins over "f"
    return out


def _pg_letters(ns):
    """p / g / f per slot, from the sign of n_pg."""
    return np.where(ns > 0, "p", np.where(ns < 0, "g", "f"))


def _pg_comp(ns):
    """Role-free p/g make-up: the three letters, p first, with every triplet
    containing an f-mode lumped into one bucket. Unlike `_pg_type` this needs no
    parent/daughter split, so it is defined for every row -- which is what
    fig_detuning needs, roles being undefined for 72% of its radial triplets."""
    s = np.sort(_pg_letters(ns), axis=1)[:, ::-1]  # sorts f < g < p, so reverse
    key = np.char.add(np.char.add(s[:, 0], s[:, 1]), s[:, 2])
    return np.where((s == "f").any(axis=1), "f", key)


# Role-free counterpart of PG_TYPE_STYLE, same colour for the same make-up.
# BLUE has no analogue: it marks mixed *parents* there, which needs roles.
PG_COMP_STYLE = {
    "ppp": ("0.55", "three $p$-modes"),
    "ppg": (GREEN, "two $p$-modes, one $g$-mode"),
    "pgg": (VERM, "one $p$-mode, two $g$-modes"),
    "ggg": ("0.05", "three $g$-modes"),
    "f": (SKY, "contains an $f$-mode ($n_{pg} = 0$)"),
}


def _pg_param_type(ns):
    """Parametric make-up, "parent|daughter pair" -- e.g. "g|gg". Only
    meaningful where slot a is the lone parent, i.e. channel == "parametric";
    MW23's five types do not apply there, they assume two parents and one
    daughter."""
    lt = _pg_letters(ns)
    pair = np.sort(lt[:, 1:], axis=1)[:, ::-1]
    return np.char.add(np.char.add(lt[:, 0], "|"),
                       np.char.add(pair[:, 0], pair[:, 1]))


def _pg_bars(ax, sub, cls, style, title):
    """Horizontal census of `cls` over `sub`. Bars count radial triplets -- the
    distinct physical triads, m being degenerate without rotation -- with the
    m-resolved count annotated. `style` maps class -> (colour, label)."""
    rad = ~sub.duplicated(["l_a", "n_a", "l_b", "n_b", "l_c", "n_c"]).to_numpy()
    # Every class in `style`, not just the populated ones: an empty class is a
    # result (no direct triplet here has two g-mode parents and a p daughter).
    names = sorted(style, key=lambda k: -int((cls == k)[rad].sum()))
    n_rad = [int((cls == k)[rad].sum()) for k in names]
    n_m = [int((cls == k).sum()) for k in names]
    y = np.arange(len(names))[::-1]
    ax.barh(y, n_rad, color=[style[k][0] for k in names], height=0.62)
    for yi, r, m in zip(y, n_rad, n_m):
        ax.text(r + max(n_rad) * 0.02, yi, f"{r}  ({m} with $m$)",
                va="center", fontsize=8, color="0.25")
    ax.set(yticks=y, yticklabels=[style[k][1] for k in names],
           xlim=(0, max(n_rad) * 1.42), xlabel="radial triplets", title=title)
    ax.grid(axis="y", visible=False)


# Compact class names for the census; the scatter figures carry the long ones.
PG_DIRECT_CENSUS = {
    "ppp": ("0.55", "$pp$ parents $\\to$ $p$"),
    "ppg": (GREEN, "$pp$ parents $\\to$ $g$"),
    "ggp": (VERM, "$gg$ parents $\\to$ $p$"),
    "ggg": ("0.05", "$gg$ parents $\\to$ $g$"),
    "mixed": (BLUE, "$pg$ parents $\\to$ either"),
    "f": (SKY, "$f$-mode in triplet"),
}


def _pg_stack(ax, x, comp, bins, **kw):
    """Stacked histogram of `x` split by `comp`, populated classes only, in
    PG_COMP_STYLE order with counts on the labels. Returns nothing."""
    order = [k for k in PG_COMP_STYLE if (comp == k).any()]
    masks = [comp == k for k in order]
    ax.hist([x[m] for m in masks], bins=bins, stacked=True,
            color=[PG_COMP_STYLE[k][0] for k in order],
            label=[f"{PG_COMP_STYLE[k][1]}  [{int(m.sum())}]"
                   for k, m in zip(order, masks)], **kw)


# gamma is looked up by frequency (coupling/modes.py DampingRates), so a mode
# with no nonadiabatic counterpart inside the tolerance gets gamma = NaN, and
# mu with it. Those rows must be dropped, not merely left undrawn: numpy sorts
# NaN last, so a descending mu ranking would hand them rank 1.
GAMMA_MU_COLS = ["mu_a", "mu_b", "mu_c", "gamma_a", "gamma_b", "gamma_c"]


def _drop_nonfinite(obs, tag=""):
    """(frame without non-finite gamma/mu rows, number dropped)."""
    ok = np.isfinite(obs[GAMMA_MU_COLS].to_numpy()).all(axis=1)
    n = int((~ok).sum())
    if n:
        print(f"  dropped {n} of {len(obs)} rows with non-finite gamma/mu"
              f"{f' [{tag}]' if tag else ''}")
        obs = obs[ok].reset_index(drop=True)
    return obs, n


def _daughter_slot(obs):
    """(daughter_index, slot to report). The slot is the daughter wherever a
    direct channel names one; elsewhere it falls back to the strongest-responding
    of b, c, which has no physical reading -- hence idx is returned alongside so
    callers can tell the two apart."""
    idx = daughter_index(obs)
    mus = obs[["mu_a", "mu_b", "mu_c"]].to_numpy()
    return idx, np.where(idx >= 0, np.maximum(idx, 0), 1 + mus[:, 1:].argmax(axis=1))


def _read_observables(usecols=None):
    """The coupling table SOURCE names, truncated at L_MAX_CUT.

    Every figure that plots triplets must come through here, or a ladder rung
    silently renders the narrow net instead of its own -- which is exactly what
    fig_eth and fig_mu used to do, reading the CSV path themselves.
    """
    path = (OUT / "lconv" / f"observables_{TAG}.csv" if SOURCE == "lconv"
            else OUT / f"observables_{TAG}.csv")
    need = None if usecols is None else sorted(set(usecols) | {"l_b", "l_c"})
    obs = pd.read_csv(path, usecols=need)
    if L_MAX_CUT is not None:
        obs = obs[obs[["l_b", "l_c"]].max(axis=1) <= L_MAX_CUT].reset_index(drop=True)
    return obs


@functools.cache
def _channel_table():
    """Ranked table plus the channel each triplet can run and the mode it
    drives, carried in the `_t` columns: the single daughter for the direct
    classes, the more strongly responding of the two for the parametric one.
    Channels with no daughter at all -- all-driven, all-damped, inactive --
    fall back to the stronger of slots b, c, which has no physical reading.

    SOURCE picks the table: "narrow" is the l <= 3 net, "lconv" the wide
    parent-restricted sweep, which is the only one carrying l > 3. L_MAX_CUT
    then truncates it on the daughter slots, which is how one sweep stands in
    for a whole ladder of nets."""
    obs, _ = _drop_nonfinite(_read_observables())
    idx, j = _daughter_slot(obs)
    mus = obs[["mu_a", "mu_b", "mu_c"]].to_numpy()
    w = obs[["omega_a", "omega_b", "omega_c"]].abs().to_numpy()
    gam = np.abs(obs[["gamma_a", "gamma_b", "gamma_c"]].to_numpy())
    ls = obs[["l_a", "l_b", "l_c"]].to_numpy()
    ns = obs[["n_a", "n_b", "n_c"]].to_numpy()
    rows = np.arange(len(obs))
    w_t, gam_t = w[rows, j], gam[rows, j]
    return obs.assign(
        channel=classify_frame(obs),
        mu_t=mus[rows, j], f_t=w_t / CD, l_t=ls[rows, j], n_t=ns[rows, j],
        pg_type=_pg_type(ns, idx, j, rows),
        abs_kappa=obs.kappa.abs(),
        delta_over_omega_t=(obs.delta.abs() / w_t),
        gamma_over_omega_t=gam_t / w_t,
        delta_over_gamma_t=obs.delta.abs() / gam_t,
        frac_detuning=(obs.delta / obs.omega_a).abs(),
    )


def _mu_panels(key, style, bg_label, title, name, sub="c"):
    """MW23 Fig. 5: mu against the four quantities it is built from, one colour
    per value of `key`, everything else in a faint background layer. `sub` is
    the subscript naming the driven mode in the axis labels. Their kappa
    abscissa is linear, which works under their mu > 1e3 cut; the full sample
    spans six decades, so it is log here."""
    df = _channel_table()
    panels = (("abs_kappa", r"$|\kappa_{abc}|$"),
              ("delta_over_omega_t", rf"$|\Delta_{{abc}}| / \omega_{sub}$"),
              ("gamma_over_omega_t", rf"$|\gamma_{sub}| / \omega_{sub}$"),
              ("delta_over_gamma_t", rf"$|\Delta_{{abc}}| / |\gamma_{sub}|$"))
    idle = ~df[key].isin(style)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.4), sharey=True)
    for ax, (col, xlabel) in zip(axes.flat, panels):
        first = col == "abs_kappa"
        ax.scatter(df[col][idle], df.mu_t[idle], s=2, color="0.85", lw=0, rasterized=True,
                   label=f"{bg_label}  [{int(idle.sum())}]" if first else None)
        for val, (color, alpha, lbl) in style.items():
            m = (df[key] == val).to_numpy()
            ax.scatter(df[col][m], df.mu_t[m], s=4, color=color, lw=0, alpha=alpha,
                       rasterized=True, label=f"{lbl}  [{int(m.sum())}]" if first else None)
        ax.set(xscale="log", yscale="log", xlabel=xlabel)
        ax.axhline(1e3, color="0.3", ls="--", lw=1.2,
                   label=r"$\mu > 10^3$ cut" if first else None)
    axes[1, 1].axvline(1.0, color="0.3", ls=":", lw=1)
    axes[0, 0].set_ylim(bottom=MU_MIN)  # sharey=True, so one call does all four
    for ax in axes[:, 0]:
        ax.set_ylabel(rf"$\mu_{sub}$")
    axes[0, 0].legend(loc="lower left", fontsize=6.5, framealpha=0.9, markerscale=3)
    # savefig(bbox="tight") pads the canvas to fit the title, so these titles
    # must stay the same width or the figures stop being overlayable.
    fig.suptitle(title, y=1.0, fontsize=11)
    save(fig, name)


def _n_channels(df):
    """Direct and parametric counts, for the two channel figures' titles."""
    return (int(df.channel.str.startswith("direct").sum()),
            int((df.channel == "parametric").sum()))


# Parent = self-excited (gamma < 0), daughter = damped and nonlinearly excited,
# as in MW23 sec 2. Slot a is only the sum mode, so it is a parent in some
# channels and the daughter in others; the subscript is "dau" because MW23's
# daughter letter c collides with our slot c, and d is the granddaughter.
CHANNEL_STYLE_ROLES = {
    "direct-diff": (VERM, 0.5, "direct, difference: two parents "
                               "$\\to$ daughter at $|\\omega_a| - \\omega_x$"),
    "direct-sum": (BLUE, 0.5, "direct, sum: two parents "
                              "$\\to$ daughter at $\\omega_b + \\omega_c$"),
    "parametric": (GREEN, 0.5, "parametric: one parent $\\to$ two daughters "
                               "(stronger one plotted)"),
}


def fig_channels():
    """Split by channel, in MW23's role names."""
    df = _channel_table()
    n_dir, n_par = _n_channels(df)
    _mu_panels(
        "channel", CHANNEL_STYLE_ROLES, "no channel, no daughter",
        r"Triplets with coupling $\mu$ coupling channel,"
        f"{len(df)} triplets with $m$: "
        f"{n_dir} direct, {n_par} parametric\n"
        r"$\mu_{\rm dau} = |\omega_{\rm dau} \kappa_{abc}| /"
        r" (\Delta_{abc}^2 + \gamma_{\rm dau}^2)^{1/2}$;"
        r"  parent $\gamma < 0$,  daughter $\gamma > 0$",
        "fig_channels", sub=r"{\rm dau}")


# MW23 Fig. 5's gray/black/green/red/blue, in Okabe-Ito. The two grays carry
# their own alpha: at 0.5 over white "0.05" renders as ~0.5 and stops being
# distinguishable from "0.55".
PG_TYPE_STYLE = {
    "ppp": ("0.55", 0.9, "three $p$-modes"),
    "ggg": ("0.05", 0.9, "three $g$-modes"),
    "ppg": (GREEN, 0.5, "two $p$-mode parents $\\to$ $g$-mode daughter"),
    "ggp": (VERM, 0.5, "two $g$-mode parents $\\to$ $p$-mode daughter"),
    "mixed": (BLUE, 0.5, "one $p$-mode and one $g$-mode parent (daughter either)"),
}


def fig_channels_pg_types():
    """Split by the p/g make-up of the triplet. MW23's five types need two
    parents and one daughter, so only the direct channels are coloured;
    everything else stays in the background layer."""
    df = _channel_table()
    n_typed = int(df.pg_type.isin(PG_TYPE_STYLE).sum())
    _mu_panels(
        "pg_type", PG_TYPE_STYLE, "no direct daughter",
        f"Coupling strength: triplet type, {n_typed} of {len(df)} "
        f"triplets with $m$ have two parents and one daughter\n"
        r"$\mu_c = |\omega_c \kappa_{abc}| / (\Delta_{abc}^2 + \gamma_c^2)^{1/2}$;"
        r"  $p$-mode is $n_{pg} > 0$,  $g$-mode is $n_{pg} < 0$",
        "fig_channels_pg_types")


def fig_detuning():
    """Fractional detuning of the radial triplets, and the cut that bounds it."""
    df = _channel_table().drop_duplicates(["l_a", "n_a", "l_b", "n_b", "l_c", "n_c"])
    x = np.log10(df.frac_detuning.to_numpy())
    order = ["direct-diff", "direct-sum", "parametric", "all-driven", "inactive", "all-damped"]
    colors = {**{k: v[0] for k, v in CHANNEL_STYLE.items()},
              "all-driven": ORANGE, "inactive": SKY, "all-damped": "0.7"}
    bins = np.linspace(x.min(), x.max(), 46)

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.2))
    ax = axes[0]
    ax.hist([x[(df.channel == c).to_numpy()] for c in order], bins=bins, stacked=True,
            color=[colors[c] for c in order],
            label=[f"{c}  [{int((df.channel == c).sum())}]" for c in order])
    ax.axvline(np.median(x), color="k", ls="--", lw=1.2)
    ax.text(np.median(x) + 0.06, ax.get_ylim()[1] * 0.92,
            f"median ${_sci(10 ** np.median(x))}$", fontsize=9)
    ax.set(xlabel=r"$\log_{10}\, |\Delta_{abc}| / \omega_a$",
           ylabel="radial triplets", title=f"(a) {len(df)} radial triplets by channel")
    ax.legend(loc="upper left", fontsize=7, framealpha=0.9)

    # The right edge is the enumeration cut, not a physical fall-off: |Delta| <
    # 0.15 sqrt(GM/R^3) is a ceiling in |Delta|/omega_a that falls as 1/omega_a.
    ax = axes[1]
    w_a = np.abs(df.omega_a.to_numpy()) / CD
    ax.scatter(w_a, df.frac_detuning, s=3, color="0.6", lw=0, rasterized=True)
    grid = np.linspace(w_a.min(), w_a.max(), 200)
    ax.plot(grid, DETUNING_CUT * bg.omega_dyn / (grid * CD), color=VERM, lw=1.6,
            label=r"enumeration cut, $|\Delta| = 0.15\sqrt{GM/R^3}$")
    ax.set(yscale="log", xlabel=r"sum-mode frequency $\omega_a$ (c/d)",
           ylabel=r"$|\Delta_{abc}| / \omega_a$",
           title="(b) the upper edge is the search cut")
    ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
    save(fig, "fig_detuning")


def fig_detuning_pg():
    """Fractional detuning by p/g make-up, parametric against direct. Classes
    are role-free: MW23's role-bearing types are undefined for most radial
    triplets here. The grey outline is the triplets running neither channel."""
    df = _channel_table().drop_duplicates(["l_a", "n_a", "l_b", "n_b", "l_c", "n_c"])
    comp = _pg_comp(df[["n_a", "n_b", "n_c"]].to_numpy())
    x = np.log10(df.frac_detuning.to_numpy())
    # Bins span every radial triplet, not each panel's subset.
    bins = np.linspace(x.min(), x.max(), 46)
    par = (df.channel == "parametric").to_numpy()
    dir_ = df.channel.str.startswith("direct").to_numpy()
    rest = ~par & ~dir_

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.2), sharey=True)
    for lbl, ax, m, tag in (("a", axes[0], par, "parametric"),
                            ("b", axes[1], dir_, "direct")):
        ax.hist(x[rest], bins=bins, histtype="step", color="0.7", lw=1.2,
                label=f"no channel  [{int(rest.sum())}]")
        _pg_stack(ax, x[m], comp[m], bins)
        med = np.median(x[m])
        ax.axvline(med, color="k", ls="--", lw=1.2)
        ax.text(med + 0.06, ax.get_ylim()[1] * 0.92,
                f"median ${_sci(10 ** med)}$", fontsize=9)
        ax.set(xlabel=r"$\log_{10}\, |\Delta_{abc}| / \omega_a$",
               title=f"({lbl}) {int(m.sum())} {tag} radial "
                     r"triplets by $p$/$g$ make-up")
        ax.legend(loc="upper left", fontsize=7, framealpha=0.9)
    axes[0].set_ylabel("radial triplets")
    save(fig, "fig_detuning_pg")


def fig_pg_census():
    """The two role-bearing channels split by p/g make-up, side by side. They
    need different classifications and that is the point: direct has two
    parents and one daughter, so MW23's five types apply; parametric has one
    parent and two daughters, so its classes are parent | daughter pair."""
    df = _channel_table()
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.0))

    par = df[df.channel == "parametric"]
    ns = par[["n_a", "n_b", "n_c"]].to_numpy()
    cls = _pg_param_type(ns)
    style = {k: (PG_COMP_STYLE[c][0], f"${k[0]}$ parent $\\to$ ${k[2:]}$")
             for k, c in zip(cls, _pg_comp(ns))}
    _pg_bars(axes[0], par, cls, style,
             f"(a) parametric: one parent $\\to$ two daughters\n"
             f"{int((~par.duplicated(['l_a','n_a','l_b','n_b','l_c','n_c'])).sum())}"
             f" radial triplets")

    dr = df[df.channel.str.startswith("direct")]
    _pg_bars(axes[1], dr, dr.pg_type.to_numpy(), PG_DIRECT_CENSUS,
             f"(b) direct: two parents $\\to$ one daughter\n"
             f"{int((~dr.duplicated(['l_a','n_a','l_b','n_b','l_c','n_c'])).sum())}"
             f" radial triplets")
    fig.tight_layout()
    save(fig, "fig_pg_census")


def fig_lowfreq():
    """Which channel can reach the low frequencies, mode by mode."""
    df = _channel_table()
    f_hi = _gamma_sets()[3][0] / CD  # bottom of the linearly driven band
    key = ["l_t", "n_t"]

    direct = df[(df.channel == "direct-diff") & (df.f_t < f_hi)]
    best = direct.groupby(key).agg(f=("f_t", "first"), mu=("mu_t", "max")).reset_index()

    par = df[(df.channel == "parametric") & (df.E_th_over_E_star > 0)]
    dau = pd.concat([
        pd.DataFrame({"l_t": par[f"l_{s}"], "n_t": par[f"n_{s}"],
                      "f": par[f"omega_{s}"].abs() / CD,
                      "E_th": par.E_th_over_E_star}) for s in "bc"])
    dau = dau[dau.f < f_hi]
    floor = dau.groupby(key).agg(f=("f", "first"), E_th=("E_th", "min")).reset_index()

    fig, axes = plt.subplots(2, 1, figsize=(6.8, 6.2), sharex=True,
                             gridspec_kw=dict(hspace=0.1))
    ax = axes[0]
    ax.scatter(best.f, best.mu * 1e-12, s=26, color=VERM, lw=0)
    ax.axhspan(2e-12, 8e-12, color=ORANGE, alpha=0.25, lw=0)
    ax.text(0.02 * f_hi, 3e-12, r"ordinary combination frequencies ($\mu \approx 4$)",
            ha="left", fontsize=8, color=VERM)
    ax.set(yscale="log", ylabel=r"$A_x = \mu A_a A_b$  (at $A_{a,b} = 10^{-6}$)",
           title=f"(a) direct, difference branch: {len(best)} damped modes reachable")

    ax = axes[1]
    ax.scatter(floor.f, floor.E_th, s=26, color=GREEN, lw=0)
    ax.axhline(1e-12, color=BLUE, ls="--", lw=1.4)
    ax.text(0.98 * f_hi, 1.4e-12, "observed parent energies", ha="right",
            fontsize=8, color=BLUE)
    ax.axhline(floor.E_th.min(), color=VERM, ls=":", lw=1.2)
    ax.text(0.98 * f_hi, floor.E_th.min() * 1.4,
            f"floor $= {_sci(floor.E_th.min())}$", ha="right", fontsize=8, color=VERM)
    n_below = int((floor.E_th < 1e-12).sum())
    ax.set(yscale="log", xlabel="mode frequency (c/d)", xlim=(0, f_hi),
           ylabel=r"$\min\, E_{\rm th} / E_\star$",
           title=f"(b) parametric: {n_below} mode(s) reachable at observed parent energies"
                 if n_below else
                 f"(b) parametric: every threshold sits "
                 f"{floor.E_th.min() / 1e-12:.0f}$\\times$ above observed")
    fig.suptitle(r"Low-frequency band ($\nu < %.2f$ c/d, below the driven band)"
                 % f_hi, y=0.98, fontsize=11)
    save(fig, "fig_lowfreq")


# ---------------------------------------------------------- four-mode net
def fig_fourmode():
    df = pd.read_csv(OUT / f"four_mode_{TAG}.csv",
                     usecols=["l_d", "E_c_over_E_th"])
    grp = df.groupby("l_d").E_c_over_E_th
    stats = pd.DataFrame({"max": grp.max(), "median": grp.median()})
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(stats.index, stats["max"], "o-", color=BLUE, label="best per $l_d$")
    ax.plot(stats.index, stats["median"], "s--", color=ORANGE, lw=1.2,
            markersize=4, label="median per $l_d$")
    ax.axhline(1.0, color=VERM, ls="--", lw=1.4)
    ax.text(16.6, 1.6, "instability threshold", color=VERM, fontsize=9)
    ax.set(yscale="log", xlabel="granddaughter degree $l_d$",
           ylabel=r"$E_c / E_{\rm th}$  at  $q_{a,b} = 10^{-6}$",
           title="Four-mode systems, full net (2.2M candidates): "
                 "all far below threshold")
    ax.legend(loc="lower left", framealpha=0.9)
    save(fig, "fig_fourmode")


# ------------------------------------------------------------------- M6
def fig_m6():
    df = pd.read_csv(OUT / f"m6_network_{TAG}.csv"
                     if (OUT / f"m6_network_{TAG}.csv").exists()
                     else OUT / "m6_network.csv")
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for name, color, ls, lbl in (("E_P", BLUE, "-", "$E_P$ (parent, pumped)"),
                                 ("E_Q", SKY, "--", "$E_Q$ (parent, pumped)"),
                                 ("E_c", ORANGE, "-", "$E_c$ (daughter, direct)"),
                                 ("E_d", VERM, "-", "$E_d$ (granddaughter, parametric)")):
        ax.plot(df.t_yr, df[name], color=color, ls=ls, lw=1.2, label=lbl)
    ax.set(yscale="log", xlabel="time (yr)", ylabel=r"$E/E_\star$",
           ylim=(1e-18, 1e-4),
           title="Four-mode system driven at $3\\times$ threshold: "
                 "bounded limit cycle")
    ax.legend(loc="lower right", framealpha=0.9, fontsize=8)
    save(fig, "fig_m6")


# ------------------------------------------------------- l convergence
@functools.cache
def _lconv_table():
    """The wide-net sweep: parent-headed triplets out to l = 25, kappa at
    m = (0, 0, 0), gamma = gamma_rad + gamma_turb.

    Lives under out/lconv/ rather than beside the narrow tables so that
    _cross_models' observables_*.csv glob does not read it as another model.
    Adds `l_d`, the larger of the two daughter degrees: a triplet enters the
    net at l_max = l_d, so this column is what the cumulative slices key on.
    """
    df = pd.read_csv(OUT / "lconv" / f"observables_{TAG}.csv")
    df = df[classify_frame(df) == "parametric"].copy()
    if L_MAX_CUT is not None:
        df = df[df[["l_b", "l_c"]].max(axis=1) <= L_MAX_CUT]
    df["l_d"] = df[["l_b", "l_c"]].max(axis=1)
    df["abs_kappa"] = df.kappa.abs()
    # Routh-Hurwitz needs a_1 > 0, i.e. the triplet must dissipate on net. Where
    # the parent outruns its daughters, eq. A7's fixed point still exists but
    # nothing settles onto it, and its energy dips below E_th -- which would be
    # nonsense read as an equilibrium.
    df["settles"] = (df.gamma_a + df.gamma_b + df.gamma_c) > 0.0
    return df


def _lconv_bands(df, col, ls, pcts=(1, 10, 50)):
    """Percentiles of `col` over every triplet with l_d <= l_max, for each
    l_max in `ls`. Cumulative, because the question is what a net truncated
    here would give -- and percentiles, because the extremum is one triplet and
    moves on luck."""
    v = df[np.isfinite(df[col]) & (df[col] > 0)]
    out = np.full((len(pcts), len(ls)), np.nan)
    for i, lm in enumerate(ls):
        sub = v[col][v.l_d <= lm]
        if len(sub) > 10:
            out[:, i] = np.percentile(sub, pcts)
    return out


def fig_lconv():
    """Does the daughter net need l = 25? Read (a) and (b) together: the
    parametric threshold saturates while mu is still climbing, and they are
    answering different questions -- how hard the parent must be driven, versus
    how strongly one daughter responds once it is."""
    df = _lconv_table()
    l_top = int(df.l_d.max())
    ls = np.arange(0, l_top + 1)
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.4), sharex=True,
                             layout="constrained")

    # (a) required parent amplitude, as a band
    q = df.assign(q_th=np.sqrt(df.E_th_over_E_star.clip(lower=0)))
    p1, p10, med = _lconv_bands(q, "q_th", ls)
    ax = axes[0, 0]
    ax.fill_between(ls, p1, med, color=BLUE, alpha=0.18, lw=0,
                    label="p1 to median of the net")
    ax.plot(ls, med, "o-", color=BLUE, label="median")
    ax.plot(ls, p10, "s--", color=BLUE, lw=1.2, markersize=3.5, label="p10")
    ax.plot(ls, p1, "-", color=BLUE, lw=1.0)
    best = float(np.sqrt(df.E_th_over_E_star[df.E_th_over_E_star > 0].min()))
    ax.axhline(best, color=VERM, ls=":", lw=1.2)
    ax.text(0.98, best * 1.25, f"best in the whole net, ${_sci(best)}$",
            transform=ax.get_yaxis_transform(), ha="right", fontsize=8, color=VERM)
    ax.axhline(A_PARENT, color="0.3", ls="--", lw=1.2)
    ax.text(0.98, A_PARENT * 1.25, r"MW23 reference parent, $10^{-6}$",
            transform=ax.get_yaxis_transform(), ha="right", fontsize=8, color="0.3")
    ax.set(yscale="log", ylabel=r"$q_{\rm th} = (E_{\rm th}/E_\star)^{1/2}$",
           title="(a) parent amplitude needed to excite a pair")
    ax.legend(loc="center left", fontsize=8, framealpha=0.9)

    # (b) mu: max and median
    ax = axes[0, 1]
    mx, md = np.full(len(ls), np.nan), np.full(len(ls), np.nan)
    for i, lm in enumerate(ls):
        sub = df.mu_max[df.l_d <= lm]
        if len(sub) > 10:
            mx[i], md[i] = sub.max(), sub.median()
    ax.plot(ls, mx, "o-", color=ORANGE, label=r"$\max \mu$")
    ax.plot(ls, md, "s--", color=ORANGE, lw=1.2, markersize=3.5, label=r"median $\mu$")
    ax.set(yscale="log", ylabel=r"$\mu$", title="(b) coupling strength")

    # (c) can anything settle? Routh-Hurwitz needs net dissipation.
    ax = axes[1, 0]
    frac = np.array([df.settles[df.l_d <= lm].mean() if (df.l_d <= lm).sum() > 10
                     else np.nan for lm in ls])
    ax.plot(ls, 100 * frac, "o-", color=GREEN,
            label="strongly damped daughters supply the dissipation")
    ax.set(ylabel=r"per cent with $\sum\gamma > 0$",
           title="(c) triplets that can reach an equilibrium at all")

    # (d) supply
    ax = axes[1, 1]
    ax.plot(ls, [(df.l_d <= l).sum() for l in ls], "o-", color=VERM,
            label="parametric triplets")
    ax.plot(ls, [(df.settles & (df.l_d <= l)).sum() for l in ls], "s--",
            color=PURPLE, lw=1.2, markersize=4, label=r"of those, net damped")
    ax.set(yscale="log", ylabel="radial triplets in the net",
           title="(d) supply of resonant triplets")

    for a in axes.flat:
        a.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        if a is not axes[0, 0]:  # (a) placed its own, clear of the captions
            a.legend(loc="best", fontsize=8, framealpha=0.9)
    for a in axes[1]:  # sharex: only the bottom row carries tick labels
        a.set_xlabel(r"net truncated at $l_{\rm max}$")
    fig.suptitle(rf"Truncating the daughter net: {len(df)} parametric triplets, "
                 rf"$\gamma_{{\rm rad}}+\gamma_{{\rm turb}}$")
    save(fig, "fig_lconv")


def fig_l_ingredients():
    """Why the threshold has an l optimum. E_th goes as
    gamma_b gamma_c (1 + Delta^2/gamma^2) / kappa^2, and the three ingredients
    pull against each other. Reference slopes are the Dziembowski scalings
    quoted in the notes, not fits."""
    df = _lconv_table()
    gam = np.concatenate([df.gamma_b.to_numpy(), df.gamma_c.to_numpy()])
    l_of_gam = np.concatenate([df.l_b.to_numpy(), df.l_c.to_numpy()])
    ls = np.arange(0, int(df.l_d.max()) + 1)

    def shell(v, l_key, how):
        return np.array([getattr(v[l_key == l], how)() if (l_key == l).sum() > 10
                         else np.nan for l in ls])

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2), layout="constrained")

    # (a) detuning. The median is useless here -- it just reports half the
    # enumeration window -- so this shows the resonant tail and draws the cut.
    ax = axes[0]
    d = pd.Series(df.delta.abs().to_numpy())
    ld = df.l_d.to_numpy()
    p1 = np.array([np.percentile(d[ld == l], 1) if (ld == l).sum() > 10 else np.nan
                   for l in ls])
    p10 = np.array([np.percentile(d[ld == l], 10) if (ld == l).sum() > 10 else np.nan
                    for l in ls])
    ax.plot(ls, p10, "o-", color=BLUE, label="p10 in shell")
    ax.plot(ls, p1, "s--", color="0.55", lw=1.1, markersize=3.5, label="p1 in shell")
    cut = DETUNING_CUT * bg.omega_dyn
    ax.axhline(cut, color=VERM, ls="--", lw=1.4)
    ax.text(0.98, cut * 0.6, r"enumeration cut $0.15\sqrt{GM/R^3}$",
            transform=ax.get_yaxis_transform(), ha="right", fontsize=8, color=VERM)
    m = np.isfinite(p10) & (ls > 0)
    if m.any():
        l0, y0 = ls[m][len(ls[m]) // 2], p10[m][len(ls[m]) // 2]
        ax.plot(ls[m], y0 * (ls[m] / l0) ** -2.0, ":", color="0.3", lw=1.2,
                label=r"$l^{-2}$")
    ax.set(xlabel=r"daughter degree $l$", yscale="log",
           ylabel=r"$|\Delta_{abc}|$  $[\mathrm{rad\,s^{-1}}]$",
           title="(a) achievable detuning")

    for ax, (y_best, y_med, lbl, ttl, slope, color, best_lbl) in zip(axes[1:], (
        (shell(pd.Series(gam), l_of_gam, "min"),
         shell(pd.Series(gam), l_of_gam, "median"),
         r"$\gamma_{\rm rad}+\gamma_{\rm turb}$  $[\mathrm{s^{-1}}]$",
         r"(b) daughter damping", 2.0, VERM, "least damped"),
        (shell(df.abs_kappa, df.l_d.to_numpy(), "max"),
         shell(df.abs_kappa, df.l_d.to_numpy(), "median"),
         r"$|\kappa_{abc}|$", r"(c) coupling strength", None, GREEN,
         "most strongly coupled"),
    )):
        ax.plot(ls, y_med, "o-", color=color, label="median in shell")
        ax.plot(ls, y_best, "s--", color="0.55", lw=1.1, markersize=3.5,
                label=best_lbl)
        if slope is not None:
            m = np.isfinite(y_med) & (ls > 0)
            if m.any():
                l0, y0 = ls[m][len(ls[m]) // 2], y_med[m][len(ls[m]) // 2]
                ax.plot(ls[m], y0 * (ls[m] / l0) ** slope, ":", color="0.3", lw=1.2,
                        label=rf"$l^{{{slope:+.0f}}}$")
        ax.set(xlabel=r"daughter degree $l$", ylabel=lbl, yscale="log", title=ttl)
        _core_gnet_boundary(ax)

    for ax in axes:
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.legend(loc="best", fontsize=8, framealpha=0.9)
    save(fig, "fig_l_ingredients")


def _core_gnet_boundary(ax, y=0.03):
    """l <= 6 carries make_inlists' `core` family (|n_pg| <= 20, so p-modes and
    low-order g-modes) on top of the g-mode net; l >= 7 is g-net only, where
    every daughter is a high-order g-mode with |n| in the hundreds. Every break
    at this l is that change of mode character, not a property of l itself."""
    ax.axvline(6.5, color="0.45", ls="-.", lw=1)
    ax.text(6.8, y, "core net ends\n($p$-modes drop out)", transform=
            ax.get_xaxis_transform(), fontsize=7, color="0.45", va="bottom")


# Daughter degree bands for the amplitude figures. Coarse on purpose: per-l
# curves at 25 values are unreadable, and the question is which range matters.
L_BANDS = ((0, 3, BLUE), (4, 6, GREEN), (7, 10, ORANGE),
           (11, 15, VERM), (16, 25, PURPLE))


def fig_amplitude_sweep():
    """How hard must a parent be driven before anything happens, and which l
    answers first. The parametric panel is empty below q ~ 4e-6: that is the
    result, not a rendering fault."""
    df = _lconv_table()
    e = df[df.E_th_over_E_star > 0]
    q = np.logspace(-6.5, -2.5, 120)

    ct = _channel_table()
    dr = ct[ct.channel.str.startswith("direct")]
    l_dr = dr[["l_b", "l_c"]].max(axis=1).to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.4), layout="constrained")

    ax = axes[0]
    for lo, hi, color in L_BANDS:
        m = (e.l_d >= lo) & (e.l_d <= hi)
        if not m.any():
            continue
        eth = e.E_th_over_E_star[m].to_numpy()
        n = [(eth < qi**2).sum() for qi in q]
        ax.plot(q, n, color=color, label=rf"${lo} \leq l_d \leq {hi}$  [{int(m.sum())}]")
    ax.axvline(A_PARENT, color="0.3", ls="--", lw=1.2)
    ax.text(A_PARENT * 1.1, 0.5, r"MW23 parent $10^{-6}$", rotation=90,
            transform=ax.get_xaxis_transform(), fontsize=8, color="0.3")
    best = float(np.sqrt(e.E_th_over_E_star.min()))
    ax.axvline(best, color=VERM, ls=":", lw=1.2)
    ax.text(best * 1.1, 0.5, f"first pair, ${_sci(best)}$", rotation=90,
            transform=ax.get_xaxis_transform(), fontsize=8, color=VERM)
    ax.set(xscale="log", yscale="log", xlabel=r"parent amplitude $q_a$",
           ylabel="daughter pairs above threshold",
           title=r"(a) parametric: nothing is excitable until "
                 r"$q_a \sim 4 \times 10^{-6}$")
    ax.legend(loc="upper left", fontsize=7, framealpha=0.9)

    ax = axes[1]
    for lo, hi, color in L_BANDS:
        m = (l_dr >= lo) & (l_dr <= hi)
        if m.sum() < 10:
            continue
        mu = dr.mu_t.to_numpy()[m]
        ax.plot(q, mu.max() * q**2, color=color,
                label=rf"${lo} \leq l_d \leq {hi}$  [{int(m.sum())}]")
        ax.plot(q, np.median(mu) * q**2, color=color, ls=":", lw=1.0)
    ax.plot(q, q, color="0.3", ls="--", lw=1.2)
    ax.text(0.98, 0.02, r"$A_x = q_a$: daughter matches its parents",
            transform=ax.transAxes, ha="right", fontsize=8, color="0.3")
    ax.set(xscale="log", yscale="log", xlabel=r"parent amplitude $q_a$",
           ylabel=r"$A_x = \mu\, A_a A_b$",
           title=r"(b) direct: solid = best in band, dotted = median")
    ax.legend(loc="upper left", fontsize=7, framealpha=0.9)
    save(fig, "fig_amplitude_sweep")


def fig_daughter_energy():
    """The two ways to make a daughter, on one scale. Direct daughters are
    driven at the beat frequency and can approach their parents' amplitude;
    parametric daughters grow to the eq. A7 fixed point, where
    E_i / E_a = (omega_i/gamma_i)/(omega_a/gamma_a), i.e. roughly the parent's
    driving rate over the daughter's damping rate. Two decades separate them."""
    ct = _channel_table()
    dr = ct[ct.channel.str.startswith("direct")]
    l_dr = dr[["l_b", "l_c"]].max(axis=1).to_numpy()
    # A_x / A_a = mu A_b, so at equal parent amplitudes the ratio is mu q.
    ratio_dr = dr.mu_t.to_numpy() * A_PARENT

    par = _lconv_table()
    par = par[par.settles]
    r_a = par.omega_a / par.gamma_a
    ratio_par = np.sqrt(np.maximum(
        pd.concat([(par[f"omega_{s}"] / par[f"gamma_{s}"]) / r_a for s in "bc"]), 0))
    l_par = np.concatenate([par.l_b.to_numpy(), par.l_c.to_numpy()])

    ls = np.arange(0, 26)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.4), sharey=True,
                             layout="constrained")
    for ax, (v, lk, ttl, color, note) in zip(axes, (
        (ratio_dr, l_dr, r"(a) direct: $A_x/A_a = \mu A_b$", VERM,
         "stops at $l_d = 12$: driven parents end at $l = 6$,\n"
         r"and the triangle rule caps $l_d \leq 2 l_{\rm parent}$"),
        (np.asarray(ratio_par), l_par,
         r"(b) parametric: $A_{\rm dau}/A_a$ at the eq.~A7 fixed point", GREEN,
         "median sits near $10^{-2}$, set by\n"
         r"$\gamma_{\rm parent}/\gamma_{\rm daughter}$"),
    )):
        v = np.asarray(v)
        ok = np.isfinite(v) & (v > 0)
        mx = np.array([v[ok & (lk == l)].max() if (ok & (lk == l)).sum() > 10
                       else np.nan for l in ls])
        md = np.array([np.median(v[ok & (lk == l)]) if (ok & (lk == l)).sum() > 10
                       else np.nan for l in ls])
        ax.plot(ls, md, "o-", color=color, label="median in shell")
        ax.plot(ls, mx, "s--", color="0.55", lw=1.1, markersize=3.5,
                label="best in shell")
        ax.axhline(1.0, color="0.3", ls="--", lw=1.2)
        ax.text(0.98, 1.3, "daughter matches parent", transform=ax.get_yaxis_transform(),
                ha="right", fontsize=8, color="0.3")
        ax.text(0.03, 0.90, note, transform=ax.transAxes, fontsize=7.5,
                va="top", color="0.35")
        ax.set(yscale="log", xlabel=r"daughter degree $l$", title=ttl)
        _core_gnet_boundary(ax)
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
    axes[0].set_ylabel(r"daughter amplitude / parent amplitude")
    fig.suptitle(rf"Two routes to a daughter, at parent amplitude $10^{{-6}}$")
    save(fig, "fig_daughter_energy")


A_PARENT = 1e-6  # MW23 Fig. 6: every parent pinned here, arbitrary but chosen
MU_CUT = 1e3     # so the best-coupled daughters land near their parents


def fig_spectrum():
    """MW23 Fig. 6 -- the artificial power spectrum. Parents pinned at
    A_a = A_b = 1e-6, daughters at A_c = mu A_a A_b, triplets with mu > 1e3.

    MW23 assume every mode is kappa-unstable and driven to that amplitude,
    taking no account of linear stability. Our gamma disagrees for most of them,
    so their set is drawn faint and the subset that has a genuinely damped
    daughter is drawn solid on top."""
    df = _channel_table()
    hi = df[df.mu_t > MU_CUT]
    real = hi.channel.str.startswith("direct").to_numpy()

    ls = hi[["l_a", "l_b", "l_c"]].to_numpy()
    ns = hi[["n_a", "n_b", "n_c"]].to_numpy()
    ws = hi[["omega_a", "omega_b", "omega_c"]].abs().to_numpy() / CD
    is_dau = (ls == hi.l_t.to_numpy()[:, None]) & (ns == hi.n_t.to_numpy()[:, None])
    # Parents all sit at the same height, so one stem per distinct (l, n_pg):
    # per-triplet stems would just paint a wall over the daughters. The faint /
    # solid split is a statement about daughters, so parents are drawn once.
    par = {(int(a), int(b)): float(w)
           for a, b, w in zip(ls[~is_dau], ns[~is_dau], ws[~is_dau])}

    fig, ax = plt.subplots(figsize=(7.8, 4.6))
    ax.vlines(list(par.values()), 1e-13, A_PARENT, color=BLUE, lw=0.9, alpha=0.7,
              rasterized=True, label=f"parents, $A = 10^{{-6}}$  [{len(par)} modes]")
    for m, alpha, lw, tag in ((~real, 0.25, 0.5, "no damped daughter"),
                              (real, 0.9, 0.7, r"damped daughter ($\gamma > 0$)")):
        sub = hi[m]
        ax.vlines(sub.f_t, 1e-13, sub.mu_t * 1e-12, color=VERM, lw=lw,
                  alpha=alpha, rasterized=True,
                  label=f"daughters, {tag}  [{len(sub)}]")

    f_f = sorted(ef.omega / CD for (l, n), ef in efs.items() if n == 0)
    for i, f in enumerate(f_f):
        ax.axvline(f, color="0.25", ls=":", lw=1.2, label=None if i else
                   f"$f$-modes, {f_f[0]:.1f}-{f_f[-1]:.1f} c/d (MW23 text: $\\approx$15)")
    w_dyn = bg.omega_dyn / CD
    ax.axvline(w_dyn, color=GREEN, ls="--", lw=1.4,
               label=f"$\\sqrt{{GM/R^3}} = {w_dyn:.2f}$ c/d (MW23 caption)")

    ax.set(yscale="log", xlim=(0, 85), ylim=(3e-13, 4e-6),
           xlabel="mode frequency (c/d)",
           ylabel=r"mode amplitude $A$",
           title=f"MW23 Fig. 6: {len(hi)} triplets with $\\mu > 10^3$;  "
                 r"$A_c = \mu A_a A_b$ at $A_a = A_b = 10^{-6}$")
    ax.legend(loc="upper right", fontsize=7, framealpha=0.9)
    save(fig, "fig_spectrum")


# MW23 Fig. 7's mass colours, in Okabe-Ito.
MASS_COLOR = {2.2: VERM, 2.0: ORANGE, 1.85: BLUE, 1.7: PURPLE}
M_SUN = 1.989e33


class _ModelTable(NamedTuple):
    tag: str
    mass: float          # Msun
    style: dict          # colour + linestyle
    obs: pd.DataFrame
    idx: np.ndarray      # daughter_index
    j: np.ndarray        # reported slot
    n_dropped: int       # rows lost to non-finite gamma/mu


def _cross_models() -> list[_ModelTable]:
    """One entry per tag with a table under OUT. Mass comes from the model's own
    GYRE summary rather than shell-scripts/config.sh. Tags of equal mass share
    MW23's colour and are told apart by linestyle -- their sec 4.1 quotes
    log g 3.9, which is Table 1's 3.93 rounded, so which of our two M = 2.0 tags
    is "the representative model" is genuinely ambiguous."""
    out = []
    seen: dict[float, int] = {}
    for path in sorted(OUT.glob("observables_*.csv")):
        tag = path.stem[len("observables_"):]
        obs, n_dropped = _drop_nonfinite(pd.read_csv(path), tag)
        idx, j = _daughter_slot(obs)
        summary = MODEL_ROOT / tag / "gyre" / "summary_nad.h5"
        mass = (float(np.ravel(load_h5(summary)["M_star"])[0]) / M_SUN
                if summary.exists() else float("nan"))
        key = round(mass, 2)
        n = seen.get(key, 0)
        seen[key] = n + 1
        style = dict(color=MASS_COLOR.get(key, "0.5"), ls=("-", "--", ":", "-.")[n % 4])
        out.append(_ModelTable(tag, mass, style, obs, idx, j, n_dropped))
    if not out:
        raise FileNotFoundError(f"no observables_*.csv in {OUT}")
    return out


def _slot_values(obs, j, cols):
    """The reported slot's value of `cols`, e.g. ("mu_a", "mu_b", "mu_c")."""
    return obs[list(cols)].to_numpy()[np.arange(len(obs)), j]


def fig_rank():
    """MW23 Fig. 7 -- triplets rank-ordered by mu, one curve per model. Panel
    (a) is their reading, every triplet; (b) keeps only those whose daughter is
    actually damped. Colour is stellar mass, as in their figure."""
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.2), sharey=True)
    for t in _cross_models():
        mu_t = _slot_values(t.obs, t.j, ("mu_a", "mu_b", "mu_c"))
        drop = f"  $-${t.n_dropped} NaN" if t.n_dropped else ""
        for ax, m in zip(axes, (np.ones(len(mu_t), bool), t.idx >= 0)):
            v = np.sort(mu_t[m])[::-1]
            ax.plot(np.arange(1, len(v) + 1), v, lw=1.4, **t.style,
                    label=f"{t.tag}   $M = {t.mass:.2f}\\,M_\\odot$   "
                          f"[{int((v > 1e4).sum())}]{drop}")

    for ax, title in zip(axes, ("(a) every triplet, as MW23",
                                r"(b) damped daughter ($\gamma > 0$) only")):
        ax.axhline(1e4, color="0.3", ls=":", lw=1.2)
        ax.set(xscale="log", yscale="log", xlabel="rank", title=title)
    axes[0].set_ylim(bottom=MU_MIN)
    axes[0].set_ylabel(r"$\mu$")
    for ax in axes:  # the counts differ per panel, so each carries its own
        ax.legend(loc="lower left", fontsize=7, framealpha=0.9,
                  title=r"[ ] = triplets with $\mu > 10^4$", title_fontsize=7)
    save(fig, "fig_rank")


def fig_daughter():
    """MW23 Fig. 8 -- mu against the daughter's frequency and radial order, all
    models at once. Top row is their reading, every triplet; bottom keeps only
    the triplets whose daughter is genuinely damped. Colour is stellar mass, as
    in fig_rank. MW23 sec 4.2 read off this figure that the mu ~ 1e3 triplets
    often carry high-order p-mode daughters (f >~ 60 c/d, n >~ 10) and sometimes
    low-frequency g-mode ones (f <~ 15 c/d, n < 0)."""
    # sharex="col" so the two readings of the same axis line up row to row.
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.4), sharey=True, sharex="col")
    for t in _cross_models():
        mu_t = _slot_values(t.obs, t.j, ("mu_a", "mu_b", "mu_c"))
        f_t = np.abs(_slot_values(t.obs, t.j, ("omega_a", "omega_b", "omega_c"))) / CD
        n_t = _slot_values(t.obs, t.j, ("n_a", "n_b", "n_c"))
        drop = f"  $-${t.n_dropped} NaN" if t.n_dropped else ""
        for r, m in enumerate((np.ones(len(mu_t), bool), t.idx >= 0)):
            for c, x in enumerate((f_t, n_t)):
                axes[r, c].scatter(
                    x[m], mu_t[m], s=2, lw=0, alpha=0.3, rasterized=True,
                    color=t.style["color"],
                    label=f"{t.tag}   $M = {t.mass:.2f}\\,M_\\odot${drop}"
                          if (r, c) == (0, 0) else None)

    rows = ("(a) every triplet, as MW23", r"(b) damped daughter ($\gamma > 0$) only")
    for r in (0, 1):
        for c, xlabel in enumerate(("daughter frequency (c/d)",
                                    r"daughter radial order $n_{pg}$")):
            ax = axes[r, c]
            for y in (1e3, 1e4):
                ax.axhline(y, color="0.3", ls=":", lw=1.0)
            if c:  # p-modes to the right of n_pg = 0, g-modes to the left
                ax.axvline(0, color="0.3", ls="--", lw=1.0)
            ax.set(yscale="log", xlabel=xlabel,
                   title=f"{rows[r]} -- by {'radial order' if c else 'frequency'}")
        axes[r, 0].set_ylabel(r"$\mu$")
    axes[0, 0].set_ylim(bottom=MU_MIN)
    axes[0, 0].legend(loc="lower right", fontsize=7, framealpha=0.9, markerscale=5)
    fig.tight_layout()
    save(fig, "fig_daughter")


def _discover_models(root: pathlib.Path):
    """Sibling model dirs load_model can read: a detail dump and a
    nonadiabatic summary. Skips polytropes and half-finished runs."""
    return sorted(d for d in root.iterdir()
                  if d.is_dir() and (d / "gyre" / "summary_nad.h5").exists()
                  and any((d / "gyre" / "detail").glob("detail.*.h5")))


ALL = {
    "propagation": fig_propagation, "gamma": fig_gamma,
    "gamma_panels": fig_gamma_panels, "kappa_cum": fig_kappa_cum,
    "xi_kappa_p": fig_xi_kappa_p, "xi_kappa_g": fig_xi_kappa_g,
    "mu": fig_mu, "eth": fig_eth, "channels": fig_channels,
    "channels_pg_types": fig_channels_pg_types,
    "detuning": fig_detuning, "lowfreq": fig_lowfreq,
    "eth_pg": fig_eth_pg, "detuning_pg": fig_detuning_pg,
    "pg_census": fig_pg_census, "spectrum": fig_spectrum,
    "fourmode": fig_fourmode, "m6": fig_m6,
    "lconv": fig_lconv, "l_ingredients": fig_l_ingredients,
    "amplitude_sweep": fig_amplitude_sweep,
    "daughter_energy": fig_daughter_energy,
}
# Cross-model: one output for the whole grid, so these run once, not per tag.
CROSS = {"rank": fig_rank, "daughter": fig_daughter}
# Everything from "kappa_cum" on needs a coupling or four-mode CSV for the tag;
# these three come from the stellar model alone.
MODEL_ONLY = ("propagation", "gamma", "gamma_panels")


def _run(model, tag, out, figs, names, source="narrow", l_max_cut=None) -> None:
    global MODEL, TAG, OUT, FIGS, bg, efs, SOURCE, L_MAX_CUT

    MODEL, TAG, OUT, FIGS = model, tag, out, figs
    SOURCE, L_MAX_CUT = source, l_max_cut
    # The tables are cached per tag, so a multi-tag run must drop them --
    # and so must a ladder run, where only L_MAX_CUT changes between passes.
    _kappa_table.cache_clear()
    _channel_table.cache_clear()
    _lconv_table.cache_clear()
    FIGS.mkdir(parents=True, exist_ok=True)

    print(f"model {MODEL}  tag {TAG}  -> {FIGS}")
    print("loading narrow model...")
    bg, efs = load_model(MODEL)

    print("figures:")
    missing = []
    for name in names:
        try:
            ALL[name]()
        except FileNotFoundError as exc:
            missing.append(f"{name}: {exc}")
    for m in missing:
        print(f"  SKIPPED {m}")
    print("done ->", FIGS)


def _run_cross(out, figs, names) -> None:
    """Cross-model figures: one output for the whole grid, written to the plots
    root rather than any tag's subdir."""
    global OUT, FIGS
    if not names:
        return
    OUT, FIGS = out, figs
    FIGS.mkdir(parents=True, exist_ok=True)
    print(f"cross-model figures -> {FIGS}")
    for name in names:
        try:
            CROSS[name]()
        except FileNotFoundError as exc:
            print(f"  SKIPPED {name}: {exc}")


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", type=pathlib.Path, default=MODEL)
    ap.add_argument("--tag", default=None,
                    help="CSV suffix in out/; defaults to the model dir name")
    ap.add_argument("--all-tags", action="store_true",
                    help="every model dir beside --model with GYRE output, one "
                         "figure set each; ignores --tag")
    ap.add_argument("--out", type=pathlib.Path, default=OUT)
    ap.add_argument("--figs", type=pathlib.Path, default=None,
                    help="output dir; defaults to out/plots/<tag>, so models "
                         "never overwrite each other. With --all-tags it is a "
                         "parent and each tag still gets its own subdir")
    ap.add_argument("--only", nargs="*", choices=sorted(ALL) + sorted(CROSS),
                    default=None,
                    help="subset of figures; default is all that have inputs")
    ap.add_argument("--source", choices=("narrow", "lconv"), default="narrow",
                    help="which coupling table to plot: the l <= 3 net, or the "
                         "wide parent-restricted sweep under out/lconv/")
    ap.add_argument("--l-max-cut", type=int, default=None,
                    help="truncate the net at this daughter degree, so one "
                         "sweep stands in for a ladder of nets")
    ap.add_argument("--model-only", action="store_true",
                    help="shorthand for the figures that need no coupling "
                         f"table: --only {' '.join(MODEL_ONLY)}")
    args = ap.parse_args()

    global MODEL_ROOT
    MODEL_ROOT = args.model.parent

    names = MODEL_ONLY if args.model_only else (args.only or list(ALL) + list(CROSS))
    per_tag = [n for n in names if n in ALL]
    cross = [n for n in names if n in CROSS]
    root = args.figs or pathlib.Path("out/plots")

    if not args.all_tags:
        tag = args.tag or args.model.name
        _run(args.model, tag, args.out, args.figs or root / tag, per_tag,
             source=args.source, l_max_cut=args.l_max_cut)
        _run_cross(args.out, root, cross)
        return 0

    models = _discover_models(args.model.parent)
    if not models:
        print(f"no readable model dirs in {args.model.parent}", file=sys.stderr)
        return 1
    print(f"{len(models)} tags: {' '.join(m.name for m in models)}\n")
    for model in models:
        _run(model, model.name, args.out, root / model.name, per_tag,
             source=args.source, l_max_cut=args.l_max_cut)
        print()
    _run_cross(args.out, root, cross)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
