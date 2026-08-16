"""GYRE eigenfunctions, renormalised to 2 omega^2 int rho |xi|^2 d3x = E_star.

Conventions:
- Angular functions are unit-normalised, int |Y_lm|^2 dOmega = 1, matching the
  T form in angular.py. The renormalisation integral therefore carries no 4pi.
- A mode is (n_pg, l, m, s) with s = +-1 and omega = s * omega_nl. Both signs
  are separate modes (T3 convention: index a runs over both, kappa written
  with no explicit conjugates).
- gamma > 0 means damping, gamma < 0 driving. GYRE gives gamma = -Im(omega);
  radiative only, with convection frozen, so nothing comes out damped inside
  the driven band. `gamma_turb` adds the turbulent channel GYRE omits, kept
  separate from Eigenfunction.gamma: MW23 find gamma_tot ~ gamma_rad, and the
  channel classification downstream reads the radiative sign.
"""

from __future__ import annotations

from dataclasses import dataclass
import functools
import pathlib
import re

import numpy as np
from scipy.interpolate import CubicSpline

from .background import G, Background, load_background
from .gyre_io import load_h5

DETAIL_RE = re.compile(r"detail\.l(\d+)\.n([+-]\d+)\.h5$")


@dataclass(frozen=True)
class Eigenfunction:
    """One (l, n_pg) solution, cgs, at the positive frequency."""

    l: int
    n_pg: int
    omega: float  # rad/s
    omega_dimless: float
    gamma: float  # rad/s, > 0 damping
    xi_r: np.ndarray  # cm
    xi_h: np.ndarray  # cm
    div_xi: np.ndarray  # dimensionless
    delta_Phi: np.ndarray  # cm^2/s^2
    ddelta_Phi_dr: np.ndarray  # cm/s^2
    bg: Background  # on this mode's own grid

    @property
    def x(self) -> np.ndarray:
        return self.bg.x

    @property
    def Lambda2(self) -> float:
        return self.l * (self.l + 1)

    def inertia_integral(self) -> float:
        """int rho (xi_r^2 + Lambda^2 xi_h^2) r^2 dr, no 4pi."""
        return _inertia(self.bg, self.l, self.xi_r, self.xi_h)


@dataclass(frozen=True)
class Mode:
    n_pg: int
    l: int
    m: int
    s: int
    ef: Eigenfunction

    @property
    def omega(self) -> float:
        return self.s * self.ef.omega

    @property
    def gamma(self) -> float:
        return self.ef.gamma

    @property
    def key(self) -> tuple[int, int, int, int]:
        return (self.n_pg, self.l, self.m, self.s)

    def __repr__(self) -> str:
        return f"Mode(n={self.n_pg:+d}, l={self.l}, m={self.m:+d}, s={self.s:+d})"


def _inertia(bg: Background, l: int, xi_r: np.ndarray, xi_h: np.ndarray) -> float:
    integrand = bg.rho * (xi_r**2 + l * (l + 1) * xi_h**2) * bg.r**2
    return float(np.trapezoid(integrand, bg.r))


def nu_turb(bg: Background, omega: float) -> np.ndarray:
    """Effective turbulent viscosity, cm^2/s, from Duguid et al. (2020) eq. 13.

    nu_FIT = u_mlt l_mlt f(w) with w = |omega| / omega_c and omega_c = u_mlt /
    l_mlt, so l_mlt = v_conv / omega_conv and alpha_MLT never appears -- MESA's
    omega_conv already carries it. f is continuous at both breakpoints. It is
    positive everywhere: Duguid's negative-viscosity branch is excluded from
    the fit, which they take as the maximum estimate of the dissipation.

    Zero outside the convection zones, NaN if the model carries no MESA
    profile.
    """
    cz = (bg.v_conv > 0) & (bg.omega_conv > 0)
    w = np.divide(abs(omega), bg.omega_conv, out=np.ones_like(bg.r), where=cz)
    f = np.where(
        w < 1e-2, 5.0,
        np.where(w <= 5.0, 0.5 * w**-0.5, 0.5 * 5**1.5 * w**-2.0),
    )
    l_mlt = np.divide(bg.v_conv, bg.omega_conv, out=np.zeros_like(bg.r), where=cz)
    # The else branch is 0 where v_conv is finite and NaN where it is not.
    return np.where(cz, bg.v_conv * l_mlt * f, 0.0 * bg.v_conv)


def _viscous_F(l: int, r: np.ndarray, xi_r: np.ndarray, xi_h: np.ndarray) -> np.ndarray:
    """Higgins & Kopal (1968) dissipation function, in the form of Lai (1994)
    eq. 8.8, for xi = xi_r Y e_r + xi_h r grad Y with int |Y_lm|^2 dOmega = 1.

    F = 2 eps_ij eps_ij - (2/3) (div xi)^2 after the angular integral, so it is
    positive semi-definite. G is Lai's div xi, built from the same splined
    dxi_r/dr as the first term rather than from ef.div_xi: the radial pieces
    cancel to (1/3)(2 eps_rr - eps_hh)^2, and mixing two estimates of the same
    derivative into that cancellation can drive F negative.

    Lai's (l+|m|)!/(l-|m|)! factor is an artefact of his Jackson-normalised
    Y_lm and does not apply here, so F carries no m dependence.
    """
    L2 = l * (l + 1)
    dxi_r = CubicSpline(r, xi_r).derivative()(r)
    dxi_h = CubicSpline(r, xi_h).derivative()(r)
    trace_h = 2 * xi_r / r - L2 * xi_h / r
    return (
        2 * dxi_r**2
        + trace_h**2
        + L2 * (dxi_h + (xi_r - xi_h) / r) ** 2
        + (l - 1) * L2 * (l + 2) * (xi_h / r) ** 2
        - (2 / 3) * (dxi_r + trace_h) ** 2
    )


def gamma_turb(ef: Eigenfunction) -> float:
    """Turbulent damping rate, rad/s, from MW23 eq. 11:

        gamma = (omega^2 / E) int dr rho r^2 nu_turb F(r),

    with E = E_star because the eigenfunctions are normalised to 2 omega^2 I =
    E_star. Always positive -- turbulent viscosity only damps.

    Cut at x <= 1: _mesa_on_gyre_grid clamps GYRE's atmosphere points to the
    outermost MESA value, which lies inside the surface convection zone, so
    keeping them would smear a nonzero viscosity over the whole atmosphere.
    r > 0 because F carries 1/r terms; that point has zero measure.
    """
    bg = ef.bg
    m = (bg.r > 0) & (bg.x <= 1.0)
    r = bg.r[m]
    F = _viscous_F(ef.l, r, ef.xi_r[m], ef.xi_h[m])
    integrand = bg.rho[m] * nu_turb(bg, ef.omega)[m] * r**2 * F
    return float(ef.omega**2 * np.trapezoid(integrand, r) / bg.E_star)


def div_xi_from_lag_rho(d: dict[str, np.ndarray]) -> np.ndarray:
    """delta rho = -rho div.xi, and GYRE reports lag_rho in units of rho."""
    return -np.real(d["lag_rho"])


def div_xi_from_euler(d: dict[str, np.ndarray], bg: Background, omega: float) -> np.ndarray:
    """Horizontal momentum equation route. Undefined for l = 0, where xi_h = 0."""
    xi_r = np.real(d["xi_r"]) * bg.R
    xi_h = np.real(d["xi_h"]) * bg.R
    dPhi = np.real(d["eul_Phi"]) * G * bg.M / bg.R
    num = bg.rho * (bg.g * xi_r - omega**2 * bg.r * xi_h + dPhi)
    return num / (bg.Gamma_1 * bg.P)


def load_eigenfunctions(
    detail_dir: pathlib.Path,
    mesa_profile: pathlib.Path | None = None,
    gammas: DampingRates | None = None,
) -> dict[tuple[int, int], Eigenfunction]:
    """GYRE refines the spatial grid per mode, so each Eigenfunction carries its
    own Background. Modes that happen to share a grid share the object."""
    out: dict[tuple[int, int], Eigenfunction] = {}
    cache: dict[bytes, Background] = {}
    for path in sorted(pathlib.Path(detail_dir).glob("detail.*.h5")):
        ef = _load_one(path, mesa_profile, gammas, cache)
        out[(ef.l, ef.n_pg)] = ef
    if not out:
        raise FileNotFoundError(f"no detail files in {detail_dir}")
    return out


def _load_one(
    path: pathlib.Path,
    mesa_profile: pathlib.Path | None,
    gammas: DampingRates | None,
    cache: dict[bytes, Background],
) -> Eigenfunction:
    d = load_h5(path)
    l, n_pg = int(d["l"]), int(d["n_pg"])
    key = np.ascontiguousarray(d["x"]).tobytes()
    if key not in cache:
        cache[key] = load_background(d, mesa_profile)
    bg = cache[key]

    omega_dimless = float(np.real(d["omega"]))
    omega = omega_dimless * bg.omega_dyn
    gm_r = G * bg.M / bg.R

    xi_r = np.real(d["xi_r"]) * bg.R
    xi_h = np.real(d["xi_h"]) * bg.R
    delta_Phi = np.real(d["eul_Phi"]) * gm_r
    ddelta_Phi_dr = np.real(d["deul_Phi"]) * gm_r / bg.R

    if "lag_rho" in d:
        div_xi = div_xi_from_lag_rho(d)
    elif l > 0:
        div_xi = div_xi_from_euler(d, bg, omega)
    else:
        raise KeyError(
            f"{path.name}: l = 0 needs lag_rho (xi_h is identically zero, so the "
            "horizontal Euler route does not exist). Re-run gyre_ad.in."
        )

    # A = sqrt(E_star / (2 omega^2 I)); every field is linear in the amplitude.
    A = np.sqrt(bg.E_star / (2 * omega**2 * _inertia(bg, l, xi_r, xi_h)))

    return Eigenfunction(
        l=l,
        n_pg=n_pg,
        omega=omega,
        omega_dimless=omega_dimless,
        gamma=gammas(l, omega) if gammas is not None else np.nan,
        xi_r=xi_r * A,
        xi_h=xi_h * A,
        div_xi=div_xi * A,
        delta_Phi=delta_Phi * A,
        ddelta_Phi_dr=ddelta_Phi_dr * A,
        bg=bg,
    )


@dataclass(frozen=True)
class DampingRates:
    """gamma = -Im(omega), rad/s, positive for damping, looked up by frequency.

    Radiative damping only. Turbulent damping is neglected.

    GYRE's nonadiabatic classifier duplicates and skips n_pg
    labels for high-order g-modes, and the effect cascades -- 7 duplicates in
    1050 modes at l <= 15, but they shift every deeper label, so 524 of the
    1057 adiabatic modes end up one label off their nonadiabatic partner,
    starting at l = 4 and worsening with l. Frequency is the sound key: where
    the labels disagree the frequencies still agree to 8e-6 relative.

    Re(omega_nad) tracks omega_ad to ~1e-5 relative through the g-mode net and
    the 8-34 c/d p-band, but the shift grows with |gamma|/omega and reaches
    ~5e-3 above 40 c/d -- a fifth of the local spacing at the top of the
    p-band, hence `tol` well above the shift there rather than three orders
    inside the spacing. A NaN therefore means one of two things: the
    nonadiabatic run missed that mode, or its shift exceeded `tol` of the local
    spacing. Above 40 c/d the second is the larger effect.

    `spacing` here is the smallest adjacent gap
    """

    omega: dict[int, np.ndarray]  # per l, sorted, rad/s
    gamma: dict[int, np.ndarray]
    tol: float = 0.3

    def __call__(self, l: int, omega: float) -> float:
        w = self.omega.get(l)
        if w is None or not len(w):
            return float("nan")
        j = int(np.searchsorted(w, omega))
        j = min(max(j - 1 if j and (j == len(w) or omega - w[j - 1] < w[j] - omega) else j, 0),
                len(w) - 1)
        gaps = []
        if j > 0:
            gaps.append(w[j] - w[j - 1])
        if j + 1 < len(w):
            gaps.append(w[j + 1] - w[j])
        spacing = min(gaps) if gaps else abs(omega)
        return float(self.gamma[l][j]) if abs(w[j] - omega) <= self.tol * spacing else float("nan")


def load_gammas(summary_nad: pathlib.Path | list[pathlib.Path], omega_dyn: float) -> DampingRates:
    paths = [summary_nad] if isinstance(summary_nad, (str, pathlib.Path)) else summary_nad
    per_l: dict[int, list[tuple[float, float]]] = {}
    for path in paths:
        d = load_h5(path)
        for l, w in zip(d["l"], d["omega"]):
            per_l.setdefault(int(l), []).append(
                (float(np.real(w) * omega_dyn), float(-np.imag(w) * omega_dyn))
            )
    omega, gamma = {}, {}
    for l, pairs in per_l.items():
        a = np.array(sorted(pairs))
        omega[l], gamma[l] = a[:, 0], a[:, 1]
    return DampingRates(omega, gamma)


def build_mode_list(
    efs: dict[tuple[int, int], Eigenfunction],
    l_max: int | None = None,
    n_range: tuple[int, int] | None = None,
    require_gamma: bool = False,
) -> list[Mode]:
    """Every (n, l, m, s); omega is m-degenerate, so eigenfunctions are shared.

    `require_gamma` drops modes the nonadiabatic run missed. Off by default --
    polytropes have no gamma at all -- but on for anything that feeds mu or
    E_th, where a NaN gamma propagates silently through the whole table.
    """
    modes: list[Mode] = []
    for (l, n_pg), ef in sorted(efs.items()):
        if l_max is not None and l > l_max:
            continue
        if n_range is not None and not (n_range[0] <= n_pg <= n_range[1]):
            continue
        if require_gamma and not np.isfinite(ef.gamma):
            continue
        for m in range(-l, l + 1):
            for s in (+1, -1):
                modes.append(Mode(n_pg=n_pg, l=l, m=m, s=s, ef=ef))
    return modes


@functools.cache
def load_model(
    model_dir: str | pathlib.Path,
    detail_dir: str = "detail",
    inlist: str = "gyre_ad.in",
    nad: str | tuple[str, ...] = "summary_nad.h5",
) -> tuple[Background, dict[tuple[int, int], Eigenfunction]]:
    """Background plus every eigenfunction for models/<name>/.

    `detail_dir`, `inlist` and `nad` select an alternative run in the same
    directory -- the wide daughter net, for instance, whose two passes write
    two nonadiabatic summaries over one shared detail directory.
    """
    root = pathlib.Path(model_dir)
    gyre, mesa = root / "gyre", root / "mesa"
    detail = gyre / detail_dir
    any_detail = next(iter(sorted(detail.glob("detail.*.h5"))), None)
    if any_detail is None:
        raise FileNotFoundError(f"no detail files in {detail}")

    profile = _paired_mesa_profile(gyre / inlist, mesa)
    bg = load_background(any_detail, profile)
    found = [gyre / name for name in ((nad,) if isinstance(nad, str) else nad)]
    found = [p for p in found if p.exists()]
    gammas = load_gammas(found, bg.omega_dyn) if found else None
    return bg, load_eigenfunctions(detail, profile, gammas)


def _paired_mesa_profile(inlist: pathlib.Path, mesa_dir: pathlib.Path) -> pathlib.Path | None:
    """The profileN.data whose .GYRE twin the inlist points at, or None when the
    inlist reads something else -- a polytrope FGONG, for instance."""
    text = inlist.read_text()
    match = re.search(r"^\s*file\s*=\s*'([^']+)'", text, re.M)
    if match is None:
        raise ValueError(f"no model file in {inlist}")
    name = pathlib.Path(match.group(1)).name
    if not name.endswith(".GYRE"):
        return None
    return mesa_dir / "LOGS" / name.replace(".GYRE", "")
