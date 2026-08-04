#!/usr/bin/env python3
"""Figures for the progress report (Reading/Papers/Study/figs/).

Every panel comes from data already on disk; the only computation is the
cumulative kappa of two representative triplets (milliseconds each).
"""

from __future__ import annotations

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
from coupling.modes import load_model

CD = 7.27220522e-5  # rad/s per cycle/day
FIGS = pathlib.Path("Reading/Papers/Study/figs")
FIGS.mkdir(parents=True, exist_ok=True)

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


print("loading narrow model (142 modes)...")
bg, efs = load_model("models/dsct_M2.0")


# ---------------------------------------------------------------- propagation
def fig_propagation():
    m = bg.x > 1e-3
    x, r = bg.x[m], bg.r[m]
    dlnP_dr = -bg.rho[m] * bg.g[m] / bg.P[m]
    dlnrho_dr = bg.dlnrho_dlnr[m] / r
    N2 = -bg.g[m] * (dlnrho_dr - dlnP_dr / bg.Gamma_1[m])
    cs2 = bg.Gamma_1[m] * bg.P[m] / bg.rho[m]
    to_cd = 86400 / (2 * np.pi)

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    N = np.sqrt(np.clip(N2, 0, None)) * to_cd
    ax.plot(x, np.where(N > 0, N, np.nan), color=BLUE, label=r"$N$")
    for l, ls in ((1, ":"), (2, "--"), (3, "-")):
        Sl = np.sqrt(l * (l + 1) * cs2) / r * to_cd
        ax.plot(x, Sl, color=ORANGE, ls=ls, lw=1.2, label=rf"$S_{{{l}}}$")
    ax.axhspan(7.93, 34.2, color=GREEN, alpha=0.12, lw=0)
    ax.text(0.55, 16, "linearly driven band", color=GREEN, fontsize=9)
    ax.axhspan(1.69, 4.92, color=VERM, alpha=0.12, lw=0)
    ax.text(0.55, 2.6, r"parametric daughter band ($\omega_c/2$)",
            color=VERM, fontsize=9)
    ax.axvline(0.0685, color="0.4", ls=":", lw=1)
    ax.text(0.075, 0.28, "conv. core\nboundary", fontsize=8, color="0.35")
    ax.set(xlabel=r"$x = r/R$", ylabel="frequency (c/d)", yscale="log",
           xlim=(0, 1.0), ylim=(0.2, 300),
           title=r"Propagation diagram, $M = 2.0\,M_\odot$ $\delta$ Sct model")
    ax.legend(loc="upper right", ncols=2, framealpha=0.9)
    save(fig, "fig_propagation")


# ------------------------------------------------------------------- damping
def fig_gamma():
    nar = load_h5("models/dsct_M2.0/gyre/summary_nad.h5")
    wide = [load_h5(f"models/dsct_M2.0/gyre/{n}")
            for n in ("summary_nad_wide.h5", "summary_nad_wide_hi.h5")]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for d in wide:
        f, g = np.real(d["freq"]), -np.imag(d["omega"])
        hi = d["l"] > 3
        ax.scatter(f[hi & (g > 0)], g[hi & (g > 0)], s=4, color="0.75",
                   lw=0, rasterized=True,
                   label=r"damped, $4 \leq l \leq 25$ (daughter net)"
                   if d is wide[0] else None)
    f, g, l = np.real(nar["freq"]), -np.imag(nar["omega"]), nar["l"]
    dam, drv = g > 0, g < 0
    ax.scatter(f[dam], g[dam], s=14, color=BLUE, lw=0, label=r"damped, $l \leq 3$")
    ax.scatter(f[drv], -g[drv], s=22, color=VERM, marker="D", lw=0,
               label=r"driven ($\gamma < 0$), $l \leq 3$")
    ax.set(xscale="log", yscale="log", xlabel="frequency (c/d)",
           ylabel=r"$|\gamma| \; / \; \sqrt{GM/R^3}$",
           title=r"Radiative damping and driving rates")
    ax.legend(loc="lower right", framealpha=0.9, fontsize=8)
    save(fig, "fig_gamma")


# -------------------------------------------------------------- kappa_cum
def fig_kappa_cum():
    kap = pd.read_csv("out/kappa_dsct_M2.0.csv")
    kap = kap[(kap.m_a == 0) & (kap.m_b == 0) & (kap.m_c == 0)]
    allpos = (kap.n_a > 0) & (kap.n_b > 0) & (kap.n_c > 0)
    allneg = (kap.n_a < 0) & (kap.n_b < 0) & (kap.n_c < 0)
    rows = [kap[allpos].loc[kap[allpos].mu.idxmax()],
            kap[allneg].loc[kap[allneg].mu.idxmax()]]

    fig, axes = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True)
    for ax, r, color, tag in zip(axes, rows, (BLUE, VERM),
                                 ("p-mode triplet", "g-mode triplet")):
        keys = [(int(r.l_a), int(r.n_a)), (int(r.l_b), int(r.n_b)),
                (int(r.l_c), int(r.n_c))]
        res = kappa_abc(*(efs[k] for k in keys), (0, 0, 0))
        ax.plot(res.r / bg.R, res.cumulative, color=color)
        ax.axhline(0, color="0.6", lw=0.7)
        ax.axvline(0.0685, color="0.4", ls=":", lw=1)
        lbl = ", ".join(f"$({l},{n:+d})$" for l, n in keys)
        ax.set_title(f"{tag}: {lbl},  $\\kappa = {res.kappa:.3g}$", fontsize=10)
        ax.set_ylabel(r"$\kappa(<r)$")
    axes[0].text(0.0685, axes[0].get_ylim()[1] * 0.75, " conv. core boundary",
                 fontsize=8, color="0.35")
    axes[1].set_xlabel(r"$x = r/R$")
    fig.suptitle(r"Cumulative coupling integral: $p$ builds at the surface, "
                 r"$g$ at the core boundary", y=0.99, fontsize=11)
    save(fig, "fig_kappa_cum")


# ----------------------------------------------------- combination amplitudes
def fig_mu():
    obs = pd.read_csv("out/observables_dsct_M2.0.csv",
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
    obs = pd.read_csv("out/observables_dsct_M2.0.csv",
                      usecols=["gamma_a", "E_th_over_E_star"])
    eth = obs.E_th_over_E_star[(obs.gamma_a < 0) & (obs.E_th_over_E_star > 0)]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.hist(np.log10(eth), bins=60, color=BLUE, alpha=0.85)
    ax.axvline(np.log10(eth.min()), color=VERM, ls="--", lw=1.4)
    ax.text(np.log10(eth.min()) + 0.15, ax.get_ylim()[1] * 0.85,
            f"min $= {eth.min():.2g}$", color=VERM, fontsize=9)
    ax.axvline(-12, color=GREEN, ls="--", lw=1.4)
    ax.text(-11.7, ax.get_ylim()[1] * 0.45,
            "observed parent\nenergies $\\sim 10^{-12}$", color=GREEN, fontsize=9)
    ax.set(xlabel=r"$\log_{10} (E_{\rm th}/E_\star)$", ylabel="triplets",
           title="Parametric thresholds, driven parents: floor sits "
                 r"$1660\times$ above observed")
    save(fig, "fig_eth")


# ---------------------------------------------------------- four-mode net
def fig_fourmode():
    df = pd.read_csv("out/four_mode_dsct_M2.0.csv",
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
    df = pd.read_csv("out/m6_network.csv")
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


print("figures:")
fig_propagation()
fig_gamma()
fig_kappa_cum()
fig_mu()
fig_eth()
fig_fourmode()
fig_m6()
print("done ->", FIGS)
