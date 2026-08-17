import numpy as np
import pytest

from coupling.network import Coupling, Network, NetworkMode, self_coupled, three_mode

# Exact resonance, no damping: omega_a + omega_b + omega_c = 0 with the sum mode
# carrying the negative sign, as a RadialTriplet does.
W = (-5.0e-4, 3.0e-4, 2.0e-4)
KAPPA = 30.0
Q0 = np.array([1e-6, 2e-7, 1.5e-7], dtype=complex)
T_END = 4.0e7  # ~ 460 d, several coupling times


@pytest.fixture(scope="module")
def conservative():
    return three_mode(W, (0.0, 0.0, 0.0), KAPPA)


@pytest.fixture(scope="module")
def run(conservative):
    return conservative.integrate(Q0, T_END)


def test_energy_conserved(conservative, run):
    """Test 11a: sum E_x is conserved when gamma = 0 and Delta = 0."""
    _, q = run
    E = conservative.energy(q).sum(axis=0)
    assert np.max(np.abs(E - E[0])) / E[0] < 1e-8


def test_manley_rowe(conservative, run):
    """Test 11b: N_a - N_b and N_a - N_c are conserved."""
    _, q = run
    N = conservative.action(q)
    for other in (1, 2):
        d = N[0] - N[other]
        assert np.max(np.abs(d - d[0])) / np.max(np.abs(N)) < 1e-8


def test_relabelling_is_a_no_op():
    """Test 11c: which mode is called `a` cannot matter."""
    order = [2, 0, 1]
    net = three_mode(W, (0.0, 0.0, 0.0), KAPPA)
    perm = Network(
        [NetworkMode("x", W[i], 0.0) for i in order],
        [Coupling(idx=(0, 1, 2), kappa=KAPPA)],
    )
    _, q = net.integrate(Q0, T_END, n_out=200)
    _, p = perm.integrate(Q0[order], T_END, n_out=200)
    assert np.max(np.abs(p - q[order])) / np.max(np.abs(q)) < 1e-8


def test_amplitude_actually_moves(conservative, run):
    """A conservation test passes trivially if nothing happens."""
    _, q = run
    E = conservative.energy(q)
    assert np.ptp(E[1]) / E[1, 0] > 0.1


def test_self_coupled_conserves_energy():
    """b = c = d: the sum slot's factor is 1 and the repeated mode's is 2."""
    net = self_coupled(omega_a=-4.0e-4, omega_d=2.0e-4, gamma_a=0.0, gamma_d=0.0, kappa=KAPPA)
    _, q = net.integrate(np.array([1e-6, 2e-7], dtype=complex), T_END)
    E = net.energy(q).sum(axis=0)
    assert np.max(np.abs(E - E[0])) / E[0] < 1e-8
    assert np.ptp(net.energy(q)[1]) / net.energy(q)[1, 0] > 0.1


def test_rotating_frame_matches_the_lab_frame(conservative):
    """The transformation is exact, so |A| and |q| must agree."""
    t, A = conservative.integrate(Q0, 2.0e5, n_out=200)
    _, q = conservative.integrate(Q0, 2.0e5, n_out=200, rotating=False)
    assert np.max(np.abs(np.abs(A) - np.abs(q))) / np.max(np.abs(q)) < 1e-8
    assert np.max(np.abs(conservative.to_lab(t, A) - q)) / np.max(np.abs(q)) < 1e-8


def test_damping_removes_energy():
    net = three_mode(W, (0.0, 1e-8, 1e-8), KAPPA)
    _, q = net.integrate(Q0, T_END)
    E = net.energy(q).sum(axis=0)
    assert E[-1] < E[0]


def test_parametric_growth_rate_matches_theory():
    """Gamma = 2 kappa sqrt(omega_b omega_c) |q_a| for a fixed, large parent."""
    from coupling.observables import parametric_growth_rate

    q_a = 1e-5
    net = three_mode(W, (0.0, 0.0, 0.0), KAPPA)
    q0 = np.array([q_a, 1e-12, 1e-12], dtype=complex)
    expected = parametric_growth_rate(KAPPA, W[1], W[2], q_a)
    t, q = net.integrate(q0, 3.0 / expected, n_out=400)
    # Fit over the interval where the daughters are still far below the parent.
    m = np.abs(q[1]) < 1e-2 * q_a
    slope = np.polyfit(t[m], np.log(np.abs(q[1][m])), 1)[0]
    assert slope == pytest.approx(expected, rel=0.05)
