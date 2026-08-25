import numpy as np
import pytest

from coupling.network import Coupling, Network, NetworkMode, self_coupled, three_mode

# Exact resonance, no damping: omega_a = omega_b + omega_c, sum mode first, as a
# RadialTriplet orders it.
W = (5.0e-4, 3.0e-4, 2.0e-4)
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
    """Test 11b: a quantum leaves the sum mode for each one entering a pair
    mode, so N_a + N_b and N_a + N_c are conserved."""
    _, q = run
    N = conservative.action(q)
    for other in (1, 2):
        d = N[0] + N[other]
        assert np.max(np.abs(d - d[0])) / np.max(np.abs(N)) < 1e-8


def test_relabelling_is_a_no_op():
    """Test 11c: which mode is called `a` cannot matter."""
    order = [2, 0, 1]
    net = three_mode(W, (0.0, 0.0, 0.0), KAPPA)
    perm = Network(
        [NetworkMode("x", W[i], 0.0) for i in order],
        [Coupling(idx=(0, 1, 2), kappa=KAPPA, sum_slot=order.index(0))],
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
    net = self_coupled(omega_a=4.0e-4, omega_d=2.0e-4, gamma_a=0.0, gamma_d=0.0, kappa=KAPPA)
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


def test_driven_parent_settles_at_the_mw25_equilibrium():
    """MW25 eq. 6 / A7. Sharpest check of the sum-slot convention: the
    equilibrium fixes the s factors, the conjugations and the detuning at once.
    Delta > gamma_b + gamma_c makes it a fixed point rather than a limit cycle.
    """
    from coupling.observables import equilibrium_energy

    w, g, k = (1.0, 0.6, 0.5), (-0.001, 0.01, 0.01), 1.0
    net = three_mode(w, g, k)
    delta = net.detuning(net.couplings[0])
    expected = equilibrium_energy(k, w[1], w[2], g[0], g[1], g[2], delta, 1.0)
    _, A = net.integrate(np.array([1e-2, 1e-5, 1e-5], dtype=complex), 3.0e4, n_out=600,
                         rtol=1e-9)
    assert net.energy(A)[0, -1] == pytest.approx(expected, rel=1e-4)


def _shared_daughter_network(gamma=(0.0, 0.0, 0.0, 0.0), kappa=KAPPA):
    """P -> Q + c, then c -> d + d. c is a pair mode in one triplet and the sum
    mode in the other, which no single sign per mode can express.
    """
    w = (5.0e-4, 3.0e-4, 2.0e-4, 1.0e-4)
    return Network(
        [NetworkMode(nm, wx, gx) for nm, wx, gx in zip("PQcd", w, gamma)],
        [Coupling(idx=(0, 1, 2), kappa=kappa, sum_slot=0),
         Coupling(idx=(2, 3, 3), kappa=kappa, sum_slot=0)],
    )


def test_shared_mode_conserves_energy_and_actions():
    """Undamped and exactly resonant. The single-triplet Manley-Rowe pairs do
    not hold when c belongs to both legs: with fluxes f1, f2 through the two,
    dN_P = -f1, dN_Q = +f1, dN_c = f1 - f2, dN_d = +2 f2, leaving N_P + N_Q and
    N_Q - N_c - N_d/2 invariant.
    """
    net = _shared_daughter_network()
    # ~20 coupling times 1/(2 omega kappa |q|); the invariants are exact, so
    # nothing is gained by integrating for longer.
    _, A = net.integrate(np.array([2e-5, 1e-5, 3e-6, 1e-6], dtype=complex), 4e7, n_out=400)
    E, N = net.energy(A), net.action(A)
    for x in (E.sum(axis=0), N[0] + N[1], N[1] - N[2] - 0.5 * N[3]):
        assert np.ptp(x) / abs(x[0]) < 1e-6
    assert np.ptp(E[3]) / E[3, 0] > 0.1


def test_parent_parametric_network_saturates_its_parents():
    """MW25 Fig. 3: two parents, a direct daughter, and one self-coupled
    parametric daughter shared by both parents. The parametric legs hang off
    the *parents*, which is what bounds them -- attaching them to the direct
    daughter instead leaves the parents with nothing to drain them.

    Their parameters, dimensionless time. Seeds start near the threshold: below
    it the parametric daughter decays, and from a 1e-7 seed it underflows past
    atol and is lost before the parents ever cross.
    """
    w, g = (1.000, 1.001, 2.0, 0.5), (-0.01, -0.01, 0.1, 0.1)
    net = Network(
        [NetworkMode(nm, wx, gx) for nm, wx, gx in zip("abcd", w, g)],
        [Coupling(idx=(2, 0, 1), kappa=1.0, sum_slot=0),   # a + b -> c, direct
         Coupling(idx=(0, 3, 3), kappa=1.0, sum_slot=0),   # a -> d + d
         Coupling(idx=(1, 3, 3), kappa=1.0, sum_slot=0)],  # b -> d + d
    )
    _, A = net.integrate(np.array([0.05, 0.05, 1e-4, 1e-4], dtype=complex),
                         3000.0, n_out=1500, rtol=1e-9)
    E = net.energy(A)
    half = len(E[0]) // 2
    e_th = 0.1 * 0.1 / (4.0 * 1.0**2 * 0.5 * 0.5)  # threshold_energy, Delta ~ 0
    for i in (0, 1):
        assert E[i].max() < 1.0, "parent ran away; the parametric leg is not draining it"
        assert E[i][half:].mean() == pytest.approx(e_th, rel=0.5)
        assert np.ptp(np.log10(E[i][half:])) > 4.0, "parent is flat, not limit-cycling"
