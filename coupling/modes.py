"""GYRE eigenfunctions, renormalised to 2 omega^2 int rho |xi|^2 d3x = E_star.

Conventions, fixed here and nowhere else:

- Angular functions are unit-normalised, int |Y_lm|^2 dOmega = 1, matching the
  T form in angular.py. The renormalisation integral therefore carries no 4pi.
- A mode is (n_pg, l, m, s) with s = +-1 and omega = s * omega_nl. Both signs
  are separate modes (T3 convention: index a runs over both, kappa written
  with no explicit conjugates).
- gamma > 0 means damping, gamma < 0 driving. GYRE gives gamma = -Im(omega);
  T1 uses the opposite sign.
"""

from __future__ import annotations

import dataclasses
import functools
import pathlib
import re

import numpy as np

from .background import G, Background, load_background
from .gyre_io import load_h5

DETAIL_RE = re.compile(r"detail\.l(\d+)\.n([+-]\d+)\.h5$")


@dataclasses.dataclass(frozen=True)
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


@dataclasses.dataclass(frozen=True)
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
    gammas: dict[tuple[int, int], float] | None = None,
) -> dict[tuple[int, int], Eigenfunction]:
    """GYRE refines the spatial grid per mode, so each Eigenfunction carries its
    own Background. Modes that happen to share a grid share the object."""
    out: dict[tuple[int, int], Eigenfunction] = {}
    cache: dict[bytes, Background] = {}
    for path in sorted(pathlib.Path(detail_dir).glob("detail.*.h5")):
        ef = _load_one(path, mesa_profile, gammas or {}, cache)
        out[(ef.l, ef.n_pg)] = ef
    if not out:
        raise FileNotFoundError(f"no detail files in {detail_dir}")
    return out


def _load_one(
    path: pathlib.Path,
    mesa_profile: pathlib.Path | None,
    gammas: dict[tuple[int, int], float],
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
        gamma=gammas.get((l, n_pg), np.nan),
        xi_r=xi_r * A,
        xi_h=xi_h * A,
        div_xi=div_xi * A,
        delta_Phi=delta_Phi * A,
        ddelta_Phi_dr=ddelta_Phi_dr * A,
        bg=bg,
    )


def load_gammas(summary_nad: pathlib.Path, omega_dyn: float) -> dict[tuple[int, int], float]:
    """gamma = -Im(omega), converted to rad/s. Positive is damping."""
    d = load_h5(summary_nad)
    return {
        (int(l), int(n)): float(-np.imag(w) * omega_dyn)
        for l, n, w in zip(d["l"], d["n_pg"], d["omega"])
    }


def build_mode_list(
    efs: dict[tuple[int, int], Eigenfunction],
    l_max: int | None = None,
    n_range: tuple[int, int] | None = None,
) -> list[Mode]:
    """Every (n, l, m, s); omega is m-degenerate, so eigenfunctions are shared."""
    modes: list[Mode] = []
    for (l, n_pg), ef in sorted(efs.items()):
        if l_max is not None and l > l_max:
            continue
        if n_range is not None and not (n_range[0] <= n_pg <= n_range[1]):
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
) -> tuple[Background, dict[tuple[int, int], Eigenfunction]]:
    """Background plus every eigenfunction for models/<name>/.

    `detail_dir` and `inlist` select an alternative run in the same directory --
    the |n| <= 100 probe, for instance.
    """
    root = pathlib.Path(model_dir)
    gyre, mesa = root / "gyre", root / "mesa"
    detail = gyre / detail_dir
    any_detail = next(iter(sorted(detail.glob("detail.*.h5"))), None)
    if any_detail is None:
        raise FileNotFoundError(f"no detail files in {detail}")

    profile = _paired_mesa_profile(gyre / inlist, mesa)
    bg = load_background(any_detail, profile)
    nad = gyre / "summary_nad.h5"
    gammas = load_gammas(nad, bg.omega_dyn) if nad.exists() else {}
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
