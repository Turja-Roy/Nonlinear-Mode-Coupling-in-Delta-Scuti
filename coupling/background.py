"""Equilibrium structure on GYRE's radial grid, in cgs."""

from __future__ import annotations

import dataclasses
import pathlib

import numpy as np

from .gyre_io import load_h5

G = 6.67430e-8
L_SUN = 3.828e33


@dataclasses.dataclass(frozen=True)
class Background:
    x: np.ndarray
    r: np.ndarray
    rho: np.ndarray
    P: np.ndarray
    Gamma_1: np.ndarray
    g: np.ndarray
    dg_dr: np.ndarray
    dlnrho_dlnr: np.ndarray
    dGamma1_dlnrho_s: np.ndarray
    M_r: np.ndarray
    T: np.ndarray  # K
    c_P: np.ndarray  # erg/g/K
    L_r: np.ndarray  # erg/s
    M: float
    R: float
    L: float
    n_atm_clamped: int

    @property
    def E_star(self) -> float:
        return G * self.M**2 / self.R

    @property
    def omega_dyn(self) -> float:
        return float(np.sqrt(G * self.M / self.R**3))

    @property
    def has_thermal_structure(self) -> bool:
        return bool(np.isfinite(self.c_P).all()) and self.L > 0.0

    @property
    def t_thermal(self) -> np.ndarray:
        """Time to radiate the heat content above r. Sets where the adiabatic
        eigenfunctions stop being trustworthy: omega t_thermal < 2 pi."""
        if not self.has_thermal_structure:
            raise ValueError("no thermal structure in this model (polytrope?)")
        u = self.c_P * self.T * self.rho * 4.0 * np.pi * self.r**2
        above = np.concatenate([[0.0], np.cumsum(0.5 * (u[1:] + u[:-1]) * np.diff(self.r))])
        return (above[-1] - above) / self.L

    def __len__(self) -> int:
        return len(self.x)


def load_background(
    detail: pathlib.Path | dict[str, np.ndarray],
    mesa_profile: pathlib.Path | None = None,
) -> Background:
    """`mesa_profile` supplies what GYRE does not carry. A polytrope has no
    thermal structure and constant Gamma_1, so it passes None: dGamma1_dlnrho_s
    is then exactly zero and T, c_P, L_r are NaN.

    Every detail file carries the full structure on that mode's own grid, so
    passing an already-loaded dict avoids re-reading and avoids interpolation.
    """
    d = load_h5(detail) if isinstance(detail, (str, pathlib.Path)) else detail
    x, R, M = d["x"], float(d["R_star"]), float(d["M_star"])
    r = x * R
    rho, G1 = d["rho"], d["Gamma_1"]

    # c_1 = (r/R)^3 (M/M_r), so g = G M x / (c_1 R^2) stays finite at x = 0.
    g = G * M * x / (d["c_1"] * R**2)
    dg_dr = np.where(
        r > 0,
        4 * np.pi * G * rho - 2 * g / np.where(r > 0, r, 1.0),
        (4 / 3) * np.pi * G * rho,
    )

    if mesa_profile is None:
        nan = np.full_like(x, np.nan)
        cols = {"dGamma1_dlnRho_s": np.zeros_like(x), "temperature": nan, "cp": nan,
                "luminosity": nan}
        n_clamped = 0
    else:
        cols, n_clamped = _mesa_on_gyre_grid(
            d["P"], mesa_profile, ("dGamma1_dlnRho_s", "temperature", "cp", "luminosity")
        )

    return Background(
        x=x,
        r=r,
        rho=rho,
        P=d["P"],
        Gamma_1=G1,
        g=g,
        dg_dr=dg_dr,
        dlnrho_dlnr=-d["V_2"] * x**2 / G1 - d["As"],
        dGamma1_dlnrho_s=cols["dGamma1_dlnRho_s"],
        M_r=d["M_r"],
        T=cols["temperature"],
        c_P=cols["cp"],
        L_r=cols["luminosity"] * L_SUN,
        M=M,
        R=R,
        L=float(d["L_star"]),
        n_atm_clamped=n_clamped,
    )


def _mesa_on_gyre_grid(
    P: np.ndarray, mesa_profile: pathlib.Path, columns: tuple[str, ...]
) -> tuple[dict[str, np.ndarray], int]:
    """Matched in log P, which is monotonic in both grids and tracks ionisation."""
    m = np.genfromtxt(mesa_profile, skip_header=5, names=True)
    order = np.argsort(m["logP"])
    lp_mesa = m["logP"][order]

    lp = np.log10(P)
    # GYRE's atmosphere points sit below MESA's surface pressure; np.interp
    # clamps them to the outermost MESA value rather than extrapolating into
    # the H ionisation zone, where dGamma1_dlnRho_s is largest.
    out = {c: np.interp(lp, lp_mesa, m[c][order]) for c in columns}
    return out, int((lp < lp_mesa[0]).sum())
