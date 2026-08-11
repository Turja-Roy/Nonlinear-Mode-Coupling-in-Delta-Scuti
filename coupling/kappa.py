"""Three-mode coupling coefficient, Eqs. (A55)-(A62).

Only omega^2 appears, so kappa does not depend on the sign assignment of a
RadialTriplet -- one value serves all three sign classes.

No eigenfunction is differentiated numerically. Every derivative was removed
analytically in favour of div.xi and of background quantities.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline, PchipInterpolator

from .angular import AngularFactors, angular_factors
from .background import Background
from .modes import Eigenfunction
from .triplets import RadialTriplet

GROUPS = ("A55", "A56", "A57", "A58", "A59", "A60", "A61", "A62")

_BG_FIELDS = ("r", "rho", "P", "Gamma_1", "g", "dg_dr", "dlnrho_dlnr", "dGamma1_dlnrho_s")
_EF_FIELDS = ("xi_r", "xi_h", "div_xi", "delta_Phi", "ddelta_Phi_dr")


@dataclass(frozen=True)
class Fields:
    """One mode's radial functions on the integration grid."""

    xi_r: np.ndarray
    xi_h: np.ndarray
    div_xi: np.ndarray
    delta_Phi: np.ndarray
    ddelta_Phi_dr: np.ndarray
    Lambda2: float
    omega2: float


@dataclass(frozen=True)
class KappaResult:
    kappa: float
    groups: dict[str, float]
    r: np.ndarray
    dkappa_dr: np.ndarray  # already divided by 2 E_star
    cumulative: np.ndarray  # kappa(< r)
    kappa_simpson: float
    refine: int = 1
    converged: bool = True  # False if the refinement ladder ran out

    @property
    def quadrature_residual(self) -> float:
        scale = max(abs(self.kappa), 1e-300)
        return abs(self.kappa - self.kappa_simpson) / scale

    @property
    def cancellation(self) -> float:
        """max |kappa(<r)| / |kappa(R)|. Above ~1e2 the integral is dominated by
        cancellation and the result is not trustworthy."""
        return float(np.max(np.abs(self.cumulative))) / max(abs(self.kappa), 1e-300)

    def outer_fraction(self, x_split: float = 0.8) -> float:
        """Share of |dkappa/dr| accumulated outside r > x_split R."""
        w = np.abs(self.dkappa_dr)
        total = np.trapezoid(w, self.r)
        if total == 0.0:
            return 0.0
        m = self.r > x_split * self.r[-1]
        return float(np.trapezoid(w[m], self.r[m]) / total)


def _cyc(i: int) -> tuple[int, int]:
    return (i + 1) % 3, (i + 2) % 3


def A55(bg, f, ang) -> np.ndarray:
    bracket = bg.Gamma_1 * (bg.Gamma_1 + 1.0) + bg.dGamma1_dlnrho_s
    return ang.T * bg.r**2 * bg.P * bracket * f[0].div_xi * f[1].div_xi * f[2].div_xi


def A56(bg, f, ang) -> np.ndarray:
    out = np.zeros_like(bg.r)
    for i in range(3):
        j, k = _cyc(i)
        out += f[j].div_xi * f[k].div_xi * (f[i].Lambda2 * f[i].xi_h - 4.0 * f[i].xi_r)
    return ang.T * bg.r * bg.P * bg.Gamma_1 * out


def A56_parts(bg, f, ang) -> tuple[np.ndarray, np.ndarray]:
    """A56 split into its Lambda^2 xi_h and its -4 xi_r pieces. For high-order
    g-modes the first must dominate, because Lambda xi_h >> xi_r. Diagnostic
    only; kappa itself never calls this."""
    hor = np.zeros_like(bg.r)
    rad = np.zeros_like(bg.r)
    for i in range(3):
        j, k = _cyc(i)
        pre = f[j].div_xi * f[k].div_xi
        hor += pre * f[i].Lambda2 * f[i].xi_h
        rad += pre * -4.0 * f[i].xi_r
    scale = ang.T * bg.r * bg.P * bg.Gamma_1
    return scale * hor, scale * rad


def A57(bg, f, ang) -> np.ndarray:
    drho_dlnr = bg.rho * bg.dlnrho_dlnr
    return ang.T * drho_dlnr * _g4(bg) * f[0].xi_r * f[1].xi_r * f[2].xi_r


def A58(bg, f, ang) -> np.ndarray:
    out = np.zeros_like(bg.r)
    for i in range(3):
        j, k = _cyc(i)
        out += f[i].div_xi * f[j].xi_r * f[k].xi_r
    return ang.T * bg.rho * bg.r * _g4(bg) * out


def A59(bg, f, ang) -> np.ndarray:
    w = f[0].omega2 * ang.G_a + f[1].omega2 * ang.G_b + f[2].omega2 * ang.G_c
    return -bg.rho * bg.r * f[0].xi_h * f[1].xi_h * f[2].xi_h * w


def A60(bg, f, ang) -> np.ndarray:
    F = (ang.F_a, ang.F_b, ang.F_c)
    out = np.zeros_like(bg.r)
    for i in range(3):
        j, k = _cyc(i)
        w = (f[i].omega2 - 3.0 * f[j].omega2 - 3.0 * f[k].omega2) * F[i] - 2.0 * (
            f[j].omega2 * F[j] + f[k].omega2 * F[k]
        )
        out += f[i].xi_r * f[j].xi_h * f[k].xi_h * w
    return -bg.rho * bg.r * out


def A61(bg, f, ang) -> np.ndarray:
    F = (ang.F_a, ang.F_b, ang.F_c)
    out = np.zeros_like(bg.r)
    for i in range(3):
        j, k = _cyc(i)
        w = f[j].omega2 * F[j] + f[k].omega2 * F[k] - 6.0 * f[i].omega2 * ang.T
        out += f[i].xi_h * f[j].xi_r * f[k].xi_r * w
    return bg.rho * bg.r * out


def A62(bg, f, ang) -> np.ndarray:
    out = np.zeros_like(bg.r)
    for i in range(3):
        j, k = _cyc(i)
        shell = bg.dlnrho_dlnr * f[j].xi_r * f[k].xi_r + bg.r * (
            f[j].xi_r * f[k].div_xi + f[k].xi_r * f[j].div_xi
        )
        out += shell * (bg.r * f[i].ddelta_Phi_dr + 2.0 * f[i].delta_Phi)
    return ang.T * bg.rho * out


def _g4(bg) -> np.ndarray:
    return 4.0 * bg.g + bg.r * bg.dg_dr


_GROUP_FUNCS = (A55, A56, A57, A58, A59, A60, A61, A62)

# Groups A55-A58 and A62 depend on the angular factors only through T.
_UNIT_T = AngularFactors(T=1.0, F_a=0.0, F_b=0.0, F_c=0.0, G_a=0.0, G_b=0.0, G_c=0.0, S=0.0)


REFINE_LADDER = (1, 2, 4, 8, 16, 32)
QUAD_TOL = 1e-4


def kappa_abc(
    a: Eigenfunction,
    b: Eigenfunction,
    c: Eigenfunction,
    ms: tuple[int, int, int],
    refine: int | None = None,
    tol: float = QUAD_TOL,
) -> KappaResult:
    """Slot order is the argument order; permuting the arguments re-evaluates
    everything from the eigenfunctions rather than reusing a cached integral.

    With `refine` unset, the grid is refined until the spline and Simpson
    quadratures agree to `tol`. Only near-null kappa needs more than one pass:
    dkappa/dr oscillates about zero, so the integral cancels even though the
    integrand itself is well conditioned.
    """
    ang = angular_factors(a.l, b.l, c.l, *ms)
    ladder = (refine,) if refine is not None else REFINE_LADDER
    for n in ladder:
        grid, fields = _on_grid((a, b, c), n)
        out = _integrate(grid, fields, ang, a.bg.E_star, n)
        if out.quadrature_residual <= tol:
            return out
    return replace(out, converged=refine is not None)


BASIS = ("T", "F_a", "F_b", "F_c", "G_a", "G_b", "G_c")


def basis_densities(bg, f) -> dict[str, np.ndarray]:
    """dkappa/dr split by which angular factor multiplies it.

    Every angular factor is a scalar, so kappa = sum_j <factor_j> I_j with the
    seven radial integrals I_j independent of m. Groups A55-A58 and A62 carry
    only T; A59 carries G; A60, A61 carry F, plus a -6 omega^2 T remainder.
    """
    P = []  # A60's radial parts: rho r xi_r^i xi_h^j xi_h^k
    S = []  # A61's radial parts: rho r xi_h^i xi_r^j xi_r^k
    for i in range(3):
        j, k = _cyc(i)
        P.append(bg.rho * bg.r * f[i].xi_r * f[j].xi_h * f[k].xi_h)
        S.append(bg.rho * bg.r * f[i].xi_h * f[j].xi_r * f[k].xi_r)

    w = [f[i].omega2 for i in range(3)]
    unit = _UNIT_T
    d_T = (
        A55(bg, f, unit)
        + A56(bg, f, unit)
        + A57(bg, f, unit)
        + A58(bg, f, unit)
        + A62(bg, f, unit)
        - 6.0 * sum(w[i] * S[i] for i in range(3))
    )

    hhh = -bg.rho * bg.r * f[0].xi_h * f[1].xi_h * f[2].xi_h
    out = {"T": d_T}
    for i, name in enumerate(("F_a", "F_b", "F_c")):
        j, k = _cyc(i)
        out[name] = -(
            (w[i] - 3.0 * w[j] - 3.0 * w[k]) * P[i] - 2.0 * w[i] * (P[j] + P[k])
        ) + w[i] * (S[j] + S[k])
    for i, name in enumerate(("G_a", "G_b", "G_c")):
        out[name] = w[i] * hhh
    return out


def _basis_pair(a, b, c, refine: int) -> tuple[dict[str, float], dict[str, float]]:
    """Spline and Simpson versions of the same seven integrals."""
    norm = 1.0 / (2.0 * a.bg.E_star)
    bg, f = _on_grid((a, b, c), refine)
    dens = basis_densities(bg, f)
    spline = {k: _quad(v, bg.r)[0] * norm for k, v in dens.items()}
    simp = {k: float(simpson(v, x=bg.r)) * norm for k, v in dens.items()}
    return spline, simp


def kappa_from_basis(basis: dict[str, float], ang: AngularFactors) -> float:
    return sum(basis[k] * getattr(ang, k) for k in BASIS)


def kappa_all_m(
    t: RadialTriplet,
    efs: dict[tuple[int, int], Eigenfunction],
    refine: int | None = None,
    tol: float = QUAD_TOL,
) -> tuple[dict[tuple[int, int, int], float], int]:
    """Every m combination of one triplet, from one set of radial integrals.

    The escalation is judged on kappa, not on the seven components: a near-null
    kappa can come from cancellation between components that are individually
    well resolved.
    """
    a, b, c = (efs[k] for k in t.keys)
    angs = {ms: angular_factors(a.l, b.l, c.l, *ms) for ms in t.m_combinations()}
    ladder = (refine,) if refine is not None else REFINE_LADDER
    for n in ladder:
        spline, simp = _basis_pair(a, b, c, n)
        out = {ms: kappa_from_basis(spline, ang) for ms, ang in angs.items()}
        ref = {ms: kappa_from_basis(simp, ang) for ms, ang in angs.items()}
        resid = max(
            abs(out[ms] - ref[ms]) / max(abs(out[ms]), 1e-300) for ms in out
        )
        if resid <= tol:
            return out, n
    return out, n


def prepare(a, b, c, ms, refine: int = 1):
    """(grid, fields, angular factors) — the inputs every group function takes."""
    grid, fields = _on_grid((a, b, c), refine)
    return grid, fields, angular_factors(a.l, b.l, c.l, *ms)


def kappa_for_triplet(
    t: RadialTriplet,
    efs: dict[tuple[int, int], Eigenfunction],
    ms: tuple[int, int, int],
    refine: int | None = None,
    tol: float = QUAD_TOL,
) -> KappaResult:
    """kappa_abc keyed by a RadialTriplet instead of three eigenfunctions."""
    a, b, c = (efs[k] for k in t.keys)
    return kappa_abc(a, b, c, ms, refine=refine, tol=tol)


def _quad(y: np.ndarray, r: np.ndarray) -> tuple[float, np.ndarray]:
    """Exact integral of the cubic spline through y."""
    anti = CubicSpline(r, y).antiderivative()
    cum = anti(r) - anti(r[0])
    return float(cum[-1]), cum


def _integrate(bg, f, ang: AngularFactors, E_star: float, refine: int = 1) -> KappaResult:
    norm = 1.0 / (2.0 * E_star)
    total = np.zeros_like(bg.r)
    groups: dict[str, float] = {}
    for name, func in zip(GROUPS, _GROUP_FUNCS):
        d = func(bg, f, ang)
        groups[name] = _quad(d, bg.r)[0] * norm
        total += d
    total *= norm
    kappa, cum = _quad(total, bg.r)
    return KappaResult(
        kappa=kappa,
        groups=groups,
        r=bg.r,
        dkappa_dr=total,
        cumulative=cum,
        kappa_simpson=float(simpson(total, x=bg.r)),
        refine=refine,
    )


def _cut(x: np.ndarray) -> int:
    """Index one past x = 1. The MESA grid continues into the atmosphere; T3's
    integral stops at the photosphere, and dGamma1_dlnrho_s is clamped above
    it anyway. A polytrope grid ends below 1, so this is a no-op there."""
    return int(np.searchsorted(x, 1.0, side="right"))


def _on_grid(efs, refine: int):
    """Background and eigenfunctions on one common integration grid.

    GYRE refines the spatial grid per mode, so the three modes of a triplet do
    not in general share one. The master grid is the union of theirs -- every
    point then has an exact background value from at least one of the three,
    and only the eigenfunctions are interpolated. A mode is smooth wherever
    GYRE declined to refine it, which is exactly where the other two add
    points, so cubic splines are enough.
    """
    shared = efs[0].bg is efs[1].bg is efs[2].bg
    if shared:
        sl = slice(0, _cut(efs[0].bg.x))
        x = efs[0].bg.x[sl]
        bg_c = {name: getattr(efs[0].bg, name)[sl] for name in _BG_FIELDS}
        ef_c = [{name: getattr(ef, name)[sl] for name in _EF_FIELDS} for ef in efs]
        if refine > 1:
            x_new = _subdivide(x, refine)
            bg_c = {k: PchipInterpolator(x, v)(x_new) for k, v in bg_c.items()}
            ef_c = [{k: CubicSpline(x, v)(x_new) for k, v in d.items()} for d in ef_c]
    else:
        xs = [ef.bg.x[: _cut(ef.bg.x)] for ef in efs]
        x = np.unique(np.concatenate(xs))
        if refine > 1:
            x = _subdivide(x, refine)
        # Concatenate the three structure samples, then interpolate once. Where
        # the sample already holds the point the interpolant reproduces it, so
        # only genuinely new points are approximated. PCHIP because background
        # quantities have kinks; the eigenfunctions do not.
        order = np.argsort(np.concatenate(xs))
        x_all = np.concatenate(xs)[order]
        keep = np.concatenate([[True], np.diff(x_all) > 0])
        bg_c = {}
        for name in _BG_FIELDS:
            v = np.concatenate([getattr(ef.bg, name)[: _cut(ef.bg.x)] for ef in efs])[order]
            bg_c[name] = PchipInterpolator(x_all[keep], v[keep])(x)
        ef_c = [
            {
                name: CubicSpline(xi, getattr(ef, name)[: _cut(ef.bg.x)])(x)
                for name in _EF_FIELDS
            }
            for ef, xi in zip(efs, xs)
        ]

    grid = _Grid(**bg_c)
    fields = tuple(
        Fields(**d, Lambda2=float(ef.Lambda2), omega2=ef.omega**2)
        for d, ef in zip(ef_c, efs)
    )
    return grid, fields


def _subdivide(x: np.ndarray, refine: int) -> np.ndarray:
    """refine-1 extra points inside every cell, endpoints preserved."""
    frac = np.arange(refine) / refine
    inner = (x[:-1, None] + np.diff(x)[:, None] * frac[None, :]).ravel()
    return np.append(inner, x[-1])


@dataclass(frozen=True)
class _Grid:
    r: np.ndarray
    rho: np.ndarray
    P: np.ndarray
    Gamma_1: np.ndarray
    g: np.ndarray
    dg_dr: np.ndarray
    dlnrho_dlnr: np.ndarray
    dGamma1_dlnrho_s: np.ndarray


def grid_convergence(
    a: Eigenfunction,
    b: Eigenfunction,
    c: Eigenfunction,
    ms: tuple[int, int, int],
    refines: tuple[int, ...] = (1, 2, 4),
) -> dict[int, float]:
    """kappa at each refinement, for checking that it settles."""
    return {n: kappa_abc(a, b, c, ms, refine=n).kappa for n in refines}


__all__ = [
    "BASIS",
    "GROUPS",
    "A56_parts",
    "Fields",
    "KappaResult",
    "basis_densities",
    "grid_convergence",
    "kappa_abc",
    "kappa_all_m",
    "kappa_for_triplet",
    "kappa_from_basis",
    "QUAD_TOL",
    "REFINE_LADDER",
    "prepare",
    *GROUPS,
]
