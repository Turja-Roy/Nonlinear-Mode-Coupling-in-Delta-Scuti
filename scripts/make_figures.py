#!/usr/bin/env python3
"""Figures for the progress report (Reading/Papers/Study/figs/).

Every panel comes from data already on disk; the only computation is the
cumulative kappa of two representative triplets (milliseconds each).
"""

from __future__ import annotations

import functools
import pathlib
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from coupling.kappa import kappa_abc
from coupling.modes import gamma_turb, load_model
from coupling.observables import classify_frame, target_index

CD = 7.27220522e-5  # rad/s per cycle/day
DETUNING_CUT = 0.15  # in units of sqrt(GM/R^3), as in run_stage345.py

# argparse defaults; main() sets the rest.
MODEL = pathlib.Path("models/dsct_M2.0")
OUT = pathlib.Path("out")
TAG = FIGS = bg = efs = None

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

    Returns (hi_l, lo_damped, lo_driven, band, turb), each an (omega, |gamma|)
    pair; `band` is the (lo, hi) omega span of the driven modes.
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
    w, g = w[lo][uniq], g[lo][uniq]
    dam, drv = g > 0, g < 0
    band = (w[drv].min(), w[drv].max())

    turb = np.array([[ef.omega, gamma_turb(ef)] for ef in efs.values()]).T
    return hi_l, (w[dam], g[dam]), (w[drv], -g[drv]), band, (turb[0], turb[1])


LBL_HI = r"damped, $4 \leq l \leq 25$ (daughter net)"
LBL_DAM = r"damped, $l \leq 3$"
LBL_DRV = r"driven ($\gamma < 0$), $l \leq 3$"
LBL_TRB = r"turbulent, $l \leq 3$"


def _gamma_axes(ax, band, top_axis=True):
    ax.axvspan(*band, color=VERM, alpha=0.07, lw=0)
    ax.set(xscale="log", yscale="log")
    if top_axis:
        sec = ax.secondary_xaxis("top", functions=(lambda x: x / CD,
                                                   lambda x: x * CD))
        sec.set_xlabel("frequency (c/d)", fontsize=9)


def fig_gamma():
    hi_l, dam, drv, band, trb = _gamma_sets()
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.scatter(*hi_l, s=4, color="0.75", lw=0, rasterized=True, label=LBL_HI)
    ax.scatter(*dam, s=14, color=BLUE, lw=0, label=LBL_DAM)
    ax.scatter(*drv, s=22, color=VERM, marker="D", lw=0, label=LBL_DRV)
    ax.scatter(*trb, s=12, color=PURPLE, marker="^", lw=0, label=LBL_TRB)
    _gamma_axes(ax, band)
    ax.set(xlabel=r"mode frequency $\omega$  $[\mathrm{rad\,s^{-1}}]$",
           ylabel=r"damping rate $|\gamma|$  $[\mathrm{s^{-1}}]$",
           title=r"Radiative and turbulent damping rates")
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
    hi_l, dam, drv, band, trb = _gamma_sets()
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
    obs = pd.read_csv(OUT / f"observables_{TAG}.csv",
                      usecols=["omega_a", "omega_b", "omega_c",
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
    obs = pd.read_csv(OUT / f"observables_{TAG}.csv",
                      usecols=["gamma_a", "E_th_over_E_star"])
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


# ------------------------------------------------------ coupling channels
CHANNEL_STYLE = {
    "direct-diff": (VERM, "difference: $\\omega_a - \\omega_x$ drives a damped daughter"),
    "direct-sum": (BLUE, "sum: two driven daughters drive the parent"),
    "parametric": (GREEN, "parametric: driven parent decays into two daughters"),
}


@functools.cache
def _channel_table():
    """Ranked table plus the channel each triplet can run and the mode it
    drives: the damped target for the direct classes, the more strongly
    responding daughter for the parametric one."""
    obs = pd.read_csv(OUT / f"observables_{TAG}.csv")
    idx = target_index(obs)
    mus = obs[["mu_a", "mu_b", "mu_c"]].to_numpy()
    w = obs[["omega_a", "omega_b", "omega_c"]].abs().to_numpy()
    gam = np.abs(obs[["gamma_a", "gamma_b", "gamma_c"]].to_numpy())
    ls = obs[["l_a", "l_b", "l_c"]].to_numpy()
    ns = obs[["n_a", "n_b", "n_c"]].to_numpy()
    j = np.where(idx >= 0, np.maximum(idx, 0), 1 + mus[:, 1:].argmax(axis=1))
    rows = np.arange(len(obs))
    w_t, gam_t = w[rows, j], gam[rows, j]
    return obs.assign(
        channel=classify_frame(obs),
        mu_t=mus[rows, j], f_t=w_t / CD, l_t=ls[rows, j], n_t=ns[rows, j],
        abs_kappa=obs.kappa.abs(),
        delta_over_omega_t=(obs.delta.abs() / w_t),
        gamma_over_omega_t=gam_t / w_t,
        delta_over_gamma_t=obs.delta.abs() / gam_t,
        frac_detuning=(obs.delta / obs.omega_a).abs(),
    )


def fig_channels():
    """MW23 Fig. 5 -- mu against the four quantities it is built from -- split
    by channel. Their mode c is the mode mu belongs to, so ours is the mode the
    channel drives. Their kappa abscissa is linear, which works under their
    mu > 1e3 cut; the full sample spans six decades, so it is log here."""
    df = _channel_table()
    panels = (("abs_kappa", r"$|\kappa_{abc}|$"),
              ("delta_over_omega_t", r"$|\Delta_{abc}| / \omega_c$"),
              ("gamma_over_omega_t", r"$|\gamma_c| / \omega_c$"),
              ("delta_over_gamma_t", r"$|\Delta_{abc}| / |\gamma_c|$"))

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.4), sharey=True)
    for ax, (col, xlabel) in zip(axes.flat, panels):
        first = col == "abs_kappa"
        idle = ~df.channel.isin(CHANNEL_STYLE)
        ax.scatter(df[col][idle], df.mu_t[idle], s=2, color="0.85", lw=0, rasterized=True,
                   label=f"no channel  [{int(idle.sum())}]" if first else None)
        for name, (color, lbl) in CHANNEL_STYLE.items():
            m = (df.channel == name).to_numpy()
            ax.scatter(df[col][m], df.mu_t[m], s=4, color=color, lw=0, alpha=0.5,
                       rasterized=True, label=f"{lbl}  [{int(m.sum())}]" if first else None)
        ax.set(xscale="log", yscale="log", xlabel=xlabel)
        ax.axhline(1e3, color="0.3", ls="--", lw=1.2,
                   label=r"$\mu > 10^3$ cut" if first else None)
    axes[1, 1].axvline(1.0, color="0.3", ls=":", lw=1)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\mu_c$")
    axes[0, 0].legend(loc="lower left", fontsize=6.5, framealpha=0.9, markerscale=3)
    n_dir = int(df.channel.str.startswith("direct").sum())
    fig.suptitle(f"MW23 Fig. 5 by coupling channel, {len(df)} triplets with $m$: "
                 f"{n_dir} direct, {int((df.channel == 'parametric').sum())} parametric\n"
                 r"$\mu_c = |\omega_c \kappa_{abc}| / (\Delta_{abc}^2 + \gamma_c^2)^{1/2}$,"
                 r" with $c$ the mode the channel drives", y=1.0, fontsize=11)
    save(fig, "fig_channels")


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
    ax.set(yscale="log", xlabel=r"parent frequency $\omega_a$ (c/d)",
           ylabel=r"$|\Delta_{abc}| / \omega_a$",
           title="(b) the upper edge is the search cut")
    ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
    save(fig, "fig_detuning")


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
    ax.set(yscale="log", xlabel="parametric daughter degree $l_d$",
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
                                 ("E_c", ORANGE, "-", "$E_c$ (direct daughter)"),
                                 ("E_d", VERM, "-", "$E_d$ (parametric daughter)")):
        ax.plot(df.t_yr, df[name], color=color, ls=ls, lw=1.2, label=lbl)
    ax.set(yscale="log", xlabel="time (yr)", ylabel=r"$E/E_\star$",
           ylim=(1e-18, 1e-4),
           title="Four-mode system driven at $3\\times$ threshold: "
                 "bounded limit cycle")
    ax.legend(loc="lower right", framealpha=0.9, fontsize=8)
    save(fig, "fig_m6")


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
    "detuning": fig_detuning, "lowfreq": fig_lowfreq,
    "fourmode": fig_fourmode, "m6": fig_m6,
}
# Everything from "kappa_cum" on needs Stage 3-5 / four-mode CSVs for the tag.
STAGE12 = ("propagation", "gamma", "gamma_panels")


def _run(model, tag, out, figs, names) -> None:
    global MODEL, TAG, OUT, FIGS, bg, efs

    MODEL, TAG, OUT, FIGS = model, tag, out, figs
    # The tables are cached per tag, so a multi-tag run must drop them.
    _kappa_table.cache_clear()
    _channel_table.cache_clear()
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
    ap.add_argument("--only", nargs="*", choices=sorted(ALL), default=None,
                    help="subset of figures; default is all that have inputs")
    ap.add_argument("--stage12", action="store_true",
                    help=f"shorthand for --only {' '.join(STAGE12)}")
    args = ap.parse_args()

    names = STAGE12 if args.stage12 else (args.only or list(ALL))
    root = args.figs or pathlib.Path("out/plots")

    if not args.all_tags:
        tag = args.tag or args.model.name
        _run(args.model, tag, args.out,
             args.figs or root / tag, names)
        return 0

    models = _discover_models(args.model.parent)
    if not models:
        print(f"no readable model dirs in {args.model.parent}", file=sys.stderr)
        return 1
    print(f"{len(models)} tags: {' '.join(m.name for m in models)}\n")
    for model in models:
        _run(model, model.name, args.out, root / model.name, names)
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
