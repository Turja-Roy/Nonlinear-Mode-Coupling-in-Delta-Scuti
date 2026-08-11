import pathlib

import numpy as np
import pytest

from coupling.kappa import kappa_abc
from coupling.observables import (
    channel,
    classify_frame,
    mu,
    nonadiabatic_fraction,
    nonadiabatic_shell,
    parametric_growth_rate,
    target_index,
    threshold_energy,
    threshold_energy_ceiling,
)

P_TRIPLET = [(2, 5), (2, 15), (2, 6)]
G_TRIPLET = [(1, -18), (1, -19), (2, -20)]

RANKED = pathlib.Path(__file__).resolve().parents[2] / "out/observables_dsct_M2.0.csv"


def test_mu_reduces_to_the_detuning_limit():
    """gamma << Delta is the regime MW23 report; mu is then gamma-independent."""
    a = mu(kappa=10.0, omega=1e-4, delta=1e-6, gamma=1e-12)
    assert a == pytest.approx(abs(1e-4 * 10.0) / 1e-6, rel=1e-9)


def test_threshold_matches_the_growth_rate_balance():
    """E_th is where Gamma equals sqrt(gamma_b gamma_c), at zero detuning."""
    k, wb, wc, gb, gc, E_star = 30.0, 3e-4, 2e-4, 1e-9, 2e-9, 5.7737e48
    E_th = threshold_energy(k, wb, wc, gb, gc, 0.0, E_star)
    q_a = np.sqrt(E_th / E_star)
    assert parametric_growth_rate(k, wb, wc, q_a) == pytest.approx(np.sqrt(gb * gc), rel=1e-12)


def test_threshold_rises_with_detuning():
    args = (30.0, 3e-4, 2e-4, 1e-9, 2e-9)
    assert threshold_energy(*args, 1e-7, 1.0) > threshold_energy(*args, 0.0, 1.0)


def test_ceiling_bounds_every_gamma_pair():
    """gamma_b gamma_c / (gamma_b + gamma_c)^2 <= 1/4, so no gamma pair exceeds it."""
    k, wb, wc, delta = 30.0, 3e-4, 2e-4, 1e-6
    ceiling = threshold_energy_ceiling(k, wb, wc, delta, 1.0)
    # Exact: E_th / ceiling = 4 gb gc / Delta^2 + 4 gb gc / (gb + gc)^2, and the
    # second term is what is bounded by 1. The first is the subleading piece the
    # asymptote drops.
    for gb, gc in [(1e-9, 1e-9), (1e-9, 5e-9), (1e-11, 3e-10), (1e-12, 1e-12)]:
        slack = 1.0 + 4 * gb * gc / delta**2
        assert threshold_energy(k, wb, wc, gb, gc, delta, 1.0) <= ceiling * slack * (1 + 1e-12)
    # Equal, tiny gammas saturate it; unequal ones fall short.
    assert threshold_energy(k, wb, wc, 1e-14, 1e-14, delta, 1.0) == pytest.approx(ceiling, rel=1e-6)
    assert threshold_energy(k, wb, wc, 1e-14, 1e-12, delta, 1.0) < 0.1 * ceiling


def test_threshold_is_zero_for_a_driven_daughter():
    assert threshold_energy(30.0, 3e-4, 2e-4, -1e-9, 2e-9, 0.0, 1.0) == 0.0


def test_channel_covers_every_gamma_sign_pattern():
    d, p = -1e-9, 1e-9  # driven, damped
    assert channel(d, p, p) == "parametric"
    assert channel(p, d, d) == "direct-sum"
    assert channel(d, d, p) == channel(d, p, d) == "direct-diff"
    assert channel(d, d, d) == "all-driven"
    assert channel(p, p, p) == "all-damped"
    assert channel(p, d, p) == channel(p, p, d) == "inactive"


@pytest.fixture(scope="module")
def ranked():
    if not RANKED.exists():
        pytest.skip(f"{RANKED.name} not built")
    import pandas as pd

    return pd.read_csv(RANKED)


def test_parametric_class_matches_the_driven_parent_filter(ranked):
    """The filter fig_eth uses, and the 6713 rows of the pipeline note."""
    cls = classify_frame(ranked)
    expect = (ranked.gamma_a < 0) & (ranked.gamma_b > 0) & (ranked.gamma_c > 0)
    assert ((cls == "parametric") == expect).all()
    assert int(expect.sum()) == 6713


def test_channel_is_constant_over_m(ranked):
    """gamma is m-degenerate, so the class belongs to the radial triplet."""
    keys = ["l_a", "n_a", "l_b", "n_b", "l_c", "n_c"]
    n = ranked.assign(channel=classify_frame(ranked)).groupby(keys).channel.nunique()
    assert (n == 1).all()


def test_direct_target_is_the_damped_mode(ranked):
    idx = target_index(ranked)
    gam = ranked[["gamma_a", "gamma_b", "gamma_c"]].to_numpy()
    direct = idx >= 0
    assert (gam[direct, idx[direct]] > 0).all()
    assert ((classify_frame(ranked).to_numpy() == "direct-sum") == (idx == 0)).all()


def test_thermal_time_decreases_outwards(bg):
    t = bg.t_thermal
    interior = bg.x <= 1.0
    assert np.all(np.diff(t[interior]) <= 0)
    assert t[0] / 3.156e7 > 1e5  # centre, > 1e5 yr
    assert t[interior][-1] / 3.156e7 < 1e-5  # photosphere, seconds


def test_nonadiabatic_shell_is_thin_and_outermost(bg):
    """Sun et al. put the adiabatic/direct divergence at r/R >~ 0.96."""
    for nu in (10.0, 30.0):
        m = nonadiabatic_shell(bg, 2 * np.pi * nu / 86400.0) & (bg.x <= 1.0)
        assert m.any()
        assert bg.x[m].min() > 0.96
        assert m[np.argmin(np.abs(bg.x - 1.0))]


def test_nonadiabatic_fraction_flags_p_modes_not_g_modes(efs, bg):
    """Test 12: p-mode triplets build kappa where the adiabatic assumption is
    weakest; g-mode triplets do not."""
    p = kappa_abc(*[efs[tuple(k)] for k in P_TRIPLET], (0, 0, 0))
    g = kappa_abc(*[efs[tuple(k)] for k in G_TRIPLET], (0, 0, 0))
    f_p = nonadiabatic_fraction(p, bg, efs[(2, 15)].omega)
    f_g = nonadiabatic_fraction(g, bg, efs[(1, -19)].omega)
    assert f_p > f_g
    assert f_g < 1e-3
