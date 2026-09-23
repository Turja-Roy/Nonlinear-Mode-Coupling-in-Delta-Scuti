/* Checks amplitude/stability/csv against the physics: Manley-Rowe on the toy
   triplets, the closed forms of the thresholds, and the pandas quirks of the
   CSV layer.

   ./compile.sh test_amplitude && build/test_amplitude   */
#include "amplitude.hpp"
#include "csv.hpp"
#include "stability.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

using amp::State;

static void close (double a, double b, double tol, const char* what) {
    if (std::abs(a - b) > tol * (1.0 + std::abs(b))) {
        std::fprintf(stderr, "FAIL %s: %.17g != %.17g\n", what, a, b);
        assert(false);
    }
}

static void same_state (const State& a, const State& b, const char* what) {
    assert(a.size() == b.size());
    for (size_t i=0 ; i<a.size() ; i++) {
        close(a[i].real(), b[i].real(), 1e-13, what);
        close(a[i].imag(), b[i].imag(), 1e-13, what);
    }
}

int main () {
    const State y = {{1.0e-6, -3.0e-7}, {2.0e-7, 5.0e-8}, {-4.0e-7, 1.5e-7}};
    const double W[3] = {2.0e-3, 1.1e-3, 0.9e-3}, G[3] = {-4e-9, 3e-9, 5e-9};
    const double KAPPA = 417.3, T = 1234.5;

    /* --- the in-place overloads odeint uses agree with the returning ones --- */
    {
        amp::Network n = amp::three_mode({W[0], W[1], W[2]}, {G[0], G[1], G[2]}, KAPPA);
        State buf(3);
        n.dq_dt(T, y, buf);        same_state(buf, n.dq_dt(T, y), "in-place dq_dt");
        n.dA_dt(T, y, buf);        same_state(buf, n.dA_dt(T, y), "in-place dA_dt");
        close(n.detuning(n.triplets()[0]), W[1] + W[2] - W[0], 1e-14, "detuning");
    }

    /* --- Manley-Rowe: the check that catches a wrong sign or factor.
       Undamped and on resonance, one triplet moves f quanta out of the sum
       mode and f into each pair mode. --- */
    {
        const double wb = 1.1e-3, wc = 0.9e-3;
        amp::Network n = amp::three_mode({wb + wc, wb, wc}, {0.0, 0.0, 0.0}, KAPPA);
        amp::Options o;
        o.n_out = 64;
        const auto sol = n.integrate({1e-6, 2e-7, 3e-7}, 2.0e7, o);
        assert(sol.t.size() == 64);

        const eig::ArrayXd N0 = n.action(sol.y.front()), E0 = n.energy(sol.y.front());
        double swing = 0.0;
        for (const State& s : sol.y) {
            const eig::ArrayXd N = n.action(s), E = n.energy(s);
            close(N[0] + N[1], N0[0] + N0[1], 1e-9, "N_a + N_b");
            close(N[0] + N[2], N0[0] + N0[2], 1e-9, "N_a + N_c");
            close(N[1] - N[2], N0[1] - N0[2], 1e-9, "N_b - N_c");
            close(E.sum(), E0.sum(), 1e-9, "energy conserved on resonance");
            swing = std::max(swing, std::abs(E[0] - E0[0]));
        }
        assert(swing > 0.1 * E0[0]);             // a static solution would pass trivially
    }
    {   // repeated mode: two quanta into d for every one out of a, which is
        // where a naive loop over the three slots double counts the d equation
        amp::Network n = amp::self_coupled(2.0e-3, 1.0e-3, 0.0, 0.0, 300.0);
        amp::Options o;
        o.n_out = 64;
        const auto sol = n.integrate({1e-6, 2e-7}, 2.0e7, o);
        const eig::ArrayXd N0 = n.action(sol.y.front());
        for (const State& s : sol.y) {
            const eig::ArrayXd N = n.action(s);
            close(2.0 * N[0] + N[1], 2.0 * N0[0] + N0[1], 1e-9, "2 N_a + N_d");
        }
    }

    /* --- damping alone decays each energy at exactly 2 gamma --- */
    {
        const double g[3] = {4e-7, 3e-7, 5e-7}, t_end = 2.0e7;
        amp::Network n = amp::three_mode({W[0], W[1], W[2]}, {g[0], g[1], g[2]}, 0.0);
        amp::Options o;
        o.n_out = 3;
        const eig::ArrayXd E = n.energy(n.integrate({1e-6, 1e-6, 1e-6}, t_end, o).y.back());
        for (int i=0 ; i<3 ; i++)
            close(E[i], 1e-12 * std::exp(-2.0 * g[i] * t_end), 1e-8, "linear decay");
    }

    /* --- A and q give the same |amplitude| --- */
    {
        amp::Network n = amp::three_mode({W[0], W[1], W[2]}, {G[0], G[1], G[2]}, KAPPA);
        amp::Options slow, full;
        slow.n_out = full.n_out = 5;
        full.slow = false;
        const auto A = n.integrate({1e-6, 1e-9, 1e-9}, 2.0e4, slow);
        const auto q = n.integrate({1e-6, 1e-9, 1e-9}, 2.0e4, full);
        for (size_t j=0 ; j<A.y.size() ; j++)
            for (int i=0 ; i<3 ; i++)
                close(std::abs(A.y[j][i]), std::abs(q.y[j][i]), 1e-7, "|A| == |q|");
        const State q_t = n.to_q(A.t.back(), A.y.back());
        for (int i=0 ; i<3 ; i++)
            close(std::abs(q_t[i]), std::abs(A.y.back()[i]), 1e-12, "to_q preserves |A|");
    }

    /* --- frozen modes hold; e_max truncates on the way up only --- */
    {
        amp::Network n = amp::three_mode({W[0], W[1], W[2]}, {G[0], G[1], G[2]}, KAPPA);
        amp::Options o;
        o.n_out = 9;
        o.frozen = {0};
        for (const State& s : n.integrate({1e-6, 1e-9, 1e-9}, 2.0e7, o).y)
            close(std::abs(s[0]), 1e-6, 1e-12, "frozen mode held");

        amp::Network all_driven =
            amp::three_mode({W[0], W[1], W[2]}, {-4e-9, -3e-9, -5e-9}, KAPPA);
        amp::Options cap;
        cap.n_out = 9;
        cap.e_max = 1e-14;
        assert(all_driven.integrate({1e-9, 1e-9, 1e-9}, 4.0e9, cap).t.back() < 4.0e9);
        // starting above the cap never fires: there is nothing to cross
        close(n.integrate({1e-6, 1e-9, 1e-9}, 2.0e7, cap).t.back(), 2.0e7, 1e-12,
              "cap does not fire from above");
    }

    /* --- stability.cpp against the closed forms --- */
    {
        const double d = 3e-7, E_star = 1e46, g = 4e-9;
        // on resonance mu is the ratio of the drive to the damping
        close(stab::mu(KAPPA, W[0], 0.0, G[1]), std::abs(W[0] * KAPPA / G[1]), 1e-14, "mu");
        close(stab::threshold_energy(KAPPA, W[1], W[2], g, g, 0.0, E_star),
              g * g / (4.0 * KAPPA * KAPPA * W[1] * W[2]) * E_star, 1e-14, "E_th");
        // E_eq is E_th with the triplet's total gamma in the detuning penalty
        close(stab::equilibrium_energy(KAPPA, W[1], W[2], G[0], g, g, d, E_star)
                  / stab::threshold_energy(KAPPA, W[1], W[2], g, g, d, E_star),
              (1.0 + d * d / ((G[0] + 2 * g) * (G[0] + 2 * g)))
                  / (1.0 + d * d / (4 * g * g)), 1e-14, "E_eq / E_th");
        // the ceiling is where E_th lands at |Delta| >> gamma with gamma_b = gamma_c
        close(stab::threshold_energy(KAPPA, W[1], W[2], g, g, 1e4 * g, E_star),
              stab::threshold_energy_ceiling(KAPPA, W[1], W[2], 1e4 * g, E_star),
              1e-7, "E_th ceiling");
        // a driven daughter has no threshold: the pair is already unstable
        close(stab::threshold_energy(KAPPA, W[1], W[2], -g, g, d), 0.0, 0.0, "E_th driven");

        /* The A7 fixed point puts every slot on E_i gamma_i / w_i = const, and
           needs the three w_i/gamma_i to share a sign. */
        const auto ee = stab::equilibrium_energies(KAPPA, {-W[0], W[1], W[2]},
                                                   {-4e-9, 3e-9, 5e-9}, d, E_star);
        const double w[3] = {-W[0], W[1], W[2]}, gg[3] = {-4e-9, 3e-9, 5e-9};
        for (int i=1 ; i<3 ; i++)
            close(ee[i] * gg[i] / w[i], ee[0] * gg[0] / w[0], 1e-12, "A7 flux ratio");
        assert(std::isnan(stab::equilibrium_energies(KAPPA, {-W[0], W[1], W[2]},
                                                     {4e-9, 3e-9, 5e-9}, d, E_star)[0]));

        // gamma < 0 is a parent, gamma > 0 a daughter, and slot a is the sum mode
        struct Case { double g[3]; stab::Channel ch; int slot; };
        const Case CASES[6] = {
            {{-1,  1,  1}, stab::Channel::parametric,  -1},
            {{ 1, -1, -1}, stab::Channel::direct_sum,   0},
            {{-1, -1,  1}, stab::Channel::direct_diff,  2},
            {{-1,  1, -1}, stab::Channel::direct_diff,  1},
            {{-1, -1, -1}, stab::Channel::all_driven,  -1},
            {{ 1,  1,  1}, stab::Channel::all_damped,  -1}};
        for (const Case& c : CASES) {
            assert(stab::channel(c.g[0], c.g[1], c.g[2]) == c.ch);
            assert(stab::daughter_slot(c.g[0], c.g[1], c.g[2]) == c.slot);
        }
        // a lone parent in a pair slot drives daughters outside this triplet
        assert(stab::channel(1, -1, 1) == stab::Channel::inactive);
        assert(std::string(stab::channel_name(stab::Channel::parametric)) == "parametric");
    }

    /* --- csv.cpp round-trips, and writes the pandas spellings --- */
    {
        const auto p = std::filesystem::temp_directory_path() / "amp_test.csv";
        {
            csv::Writer w(p, {"x", "flag", "name"});
            w.row(1.5, true, std::string("a"));
            w.row(std::nan(""), false, std::string("b"));
            w.row(-2.0, true, std::string("c"));
        }
        const csv::Table t = csv::read(p);
        assert(t.rows() == 3);
        assert(t.has("flag") && !t.has("nope"));
        assert(t.text("flag")[0] == "True" && t.text("name")[2] == "c");
        close(t.numbers("x")[0], 1.5, 0.0, "column x");
        assert(std::isnan(t.numbers("x")[1]));     // NaN is written as an empty field
        close(t.numbers("x")[2], -2.0, 0.0, "column x");
        close(csv::median({3.0, 1.0, std::nan(""), 2.0}), 2.0, 0.0, "median drops NaN");
        close(csv::quantile({1.0, 2.0, 3.0, 4.0}, 0.1), 1.3, 1e-15, "quantile");
        std::filesystem::remove(p);
    }

    std::printf("amplitude/stability/csv: all checks passed\n");
    return 0;
}
