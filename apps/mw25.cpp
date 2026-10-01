/* Mourabit & Weinberg (2025), ApJ 986:33, Figures 1-6, integrated in C++.

       ./mw25 --out out/mw25_cpp
       python3 scripts/mw25_plot.py --data out/mw25_cpp

   This writes one CSV per panel (t and the mode energies) plus a single
   lines.csv carrying the horizontal thresholds and equilibria. Plotting stays
   in Python; the point of running it here is that the trajectories then come
   out of the same Network the pipeline uses.

   Sign conventions. MW25 write signed frequencies whose sum is the detuning:
   a self-excited (parent) mode carries omega < 0, a damped (daughter) mode
   omega > 0 -- amp::Network now uses that convention directly, with no
   sum-slot role attached to a triplet's position. Values below are quoted
   from mw25_simplified/mw25_simplified.cpp, which the professor has already
   checked against the paper; this file differs only in going through the
   shared amp::Network instead of its own standalone integrator. */

#include "amplitude.hpp"
#include "csv.hpp"
#include "stability.hpp"

#include <CLI/CLI.hpp>

#include <cmath>
#include <cstdio>

using amp::Network;
using amp::State;
using amp::Triplet;

namespace {

constexpr double DAY = 86400.0;

std::filesystem::path OUT;
std::vector<std::array<std::string, 3>> LINES;   // panel, label, value

void line(const std::string& panel, const std::string& label, double v) {
    LINES.push_back({panel, label, csv::fmt(v)});
}

// One panel: integrate and write t plus one energy column per mode.
void panel(const std::string& name, const Network& net, State q0, double t_end,
           int n_out, double rtol, double e_max = 0.0,
           const std::vector<int>& frozen = {}, double t_scale = 1.0) {
    amp::Options opt;
    opt.n_out = n_out;
    opt.rtol = rtol;
    opt.e_max = e_max;
    opt.frozen = frozen;
    const amp::Solution sol = net.integrate(std::move(q0), t_end, opt);
    const auto& t = sol.t;
    std::vector<std::string> hdr = {"t"};
    for (const auto& m : net.modes()) hdr.push_back("E_" + m.name);

    csv::Writer w(OUT / (name + ".csv"), hdr);
    for (size_t j = 0; j < t.size(); ++j) {
        std::string row = csv::fmt(t[j] / t_scale);
        const eig::ArrayXd E = net.energy(sol.y[j]);
        for (eig::Index i = 0; i < E.size(); ++i) row += "," + csv::fmt(E[i]);
        w.row(row);
    }
    std::printf("  %-8s %5zu rows, t to %.4g%s\n", name.c_str(), t.size(),
                t.back() / t_scale, t.back() < t_end ? "  (stopped early: e_max)" : "");
}

/* ------------------------------------------------------------- Figures 1-2 */

void fig1() {
    // (a) self-coupled parent a = b at |w| = 0.5 driving daughter c at |w| = 1.0
    {
        const double omega_a = -1.0, gamma_a = -0.01;
        const double omega_c = 2.0, gamma_c = 0.1;
        Network net({{"c", omega_c, gamma_c}, {"a", omega_a, gamma_a}},
                    {Triplet{0, 1, 1, 1.0}});
        panel("fig1a", net, State(2, 1e-3), 6000.0, 8000, 1e-11);
    }
    // (b) two distinct parents: the extra degree of freedom breaks it
    {
        const double omega_a = -1.000, gamma_a = -0.010;
        const double omega_b = -1.004, gamma_b = -0.010;
        const double omega_c = 2.000, gamma_c = 0.1;
        Network net({{"c", omega_c, gamma_c}, {"a", omega_a, gamma_a}, {"b", omega_b, gamma_b}},
                    {Triplet{0, 1, 2, 1.0}});
        panel("fig1b", net, State(3, 1e-3), 6000.0, 8000, 1e-11, 1e2);
    }
}

void fig2() {
    /* For the randomization within 1% range */
    unsigned long s = 12345;
    auto jitter = [&s] {
        s = (1103515245UL * s + 12345UL) & 0x7fffffffUL;
        return 0.01 * (2.0 * (double(s) / 2147483647.0) - 1.0);
    };
    for (int N : {2, 10, 50}) {
        std::vector<amp::Mode> modes = {{"a", -0.5, -0.001}, {"b", -0.5, -0.0011}};
        std::vector<Triplet> trs;
        for (int i=0 ; i<N ; i++) {
            const double wd = 1.0 * (1 + jitter()), gd = 0.01 * (1 + jitter());
            const double kd = 1.0 * (1 + jitter());

            modes.push_back({"d" + std::to_string(i), wd, gd});
            trs.push_back(Triplet{2 + i, 0, 1, kd});
        }
        panel("fig2_N" + std::to_string(N), Network(modes, trs),
              State(size_t(N + 2), 1e-6), 11000.0, 3000, 1e-10, 1e0);
    }
}

/* ----------------------------------------------------------- Figures 3-5 */

// MW25 Section 5: parents a, b; direct daughter c; self-coupled parametric
// daughter d. |omega| and gamma in order a, b, c, d.
constexpr double OMEGA[4] = {-1.000, -1.001, 2.0, 0.5};
constexpr double GAMMA[4] = {-0.01, -0.01, 0.1, 0.1};

Network mixed(double k_direct, double k_param, const double g[4]) {
    return Network({{"a", OMEGA[0], g[0]}, {"b", OMEGA[1], g[1]},
                    {"c", OMEGA[2], g[2]}, {"d", OMEGA[3], g[3]}},
                   {Triplet{2, 0, 1, k_direct},      // a + b -> c
                    Triplet{0, 3, 3, k_param},       // a -> d + d
                    Triplet{1, 3, 3, k_param}});     // b -> d + d
}

// E_a_th/E_a_eq: MW25 Eq. (6)/A7 for parent a, pumped by its self-coupled
// daughter pair (d, d) -- the parent's threshold/equilibrium, not d's own.
void mixed_lines(const std::string& tag, const Network& net, double k_param,
                 const double g[4]) {
    const double dp = net.detuning(net.triplets()[1]);
    line(tag, "E_a_th", stab::threshold_energy(k_param, OMEGA[3], OMEGA[3], g[3], g[3], dp));
    line(tag, "E_a_eq",
         stab::equilibrium_energy(k_param, OMEGA[3], OMEGA[3], g[0], g[3], g[3], dp));
}

void fig3() {
    const Network net = mixed(1.0, 1.0, GAMMA);
    mixed_lines("3", net, 1.0, GAMMA);
    panel("fig3", net, State(4, 1e-3), 2000.0, 20000, 1e-11);
}

void fig4() {
    // (a) parametric coupling 100x stronger: E_th scales as kappa^-2
    {
        const Network net = mixed(0.5, 100.0, GAMMA);
        mixed_lines("4a", net, 100.0, GAMMA);
        panel("fig4a", net, State{3e-7, 3e-7, 3e-10, 3e-7}, 5000.0, 20000, 1e-11);
    }
    // (b) every linear rate 100x smaller: the cycle stretches out
    {
        double g[4];
        for (int i=0 ; i<4 ; i++) g[i] = 0.01 * GAMMA[i];
        const Network net = mixed(0.5, 1.0, g);
        mixed_lines("4b", net, 1.0, g);
        panel("fig4b", net, State{3e-7, 3e-7, 3e-10, 3e-7}, 500000.0, 20000, 1e-11);
    }
}

void fig5() {
    // Does q_c = mu q_a q_b survive the parametric leg? mu belongs to the
    // daughter of the direct triplet, MW25 Eq. (9).
    for (auto [tag, gpar, t_end, n_out] : {std::tuple{"fig5a", -0.01, 2000.0, 20000},
                                    std::tuple{"fig5b", -0.003, 40000.0, 95000}}) {
        double g[4] = {gpar, gpar, GAMMA[2], GAMMA[3]};
        const Network net = mixed(1.0, 1.0, g);
        line(std::string(tag).substr(3), "mu",
             stab::mu(1.0, OMEGA[2], net.detuning(net.triplets()[0]), g[2]));
        panel(tag, net, State(4, 1e-3), t_end, n_out, 1e-11);
    }
}

/* MW25 Section 6, their 2.0 Msun, Teff = 7202 K, log g = 3.80 model.
   Parents a, b and direct daughter c: l = {3, 1, 2}, n = {0, 3, 8}.
   Parametric daughter pair: l = {14, 15}, n = {-282, -12}.

   Three values have been changed from the paper's quoted numbers:
     Delta_direct = 7.6e-6 (quoted 7.6e-5)
     gamma_c = 2.8e-6 (quoted 2.8e-5)
     gamma_d = {1.4, 1.5}e-5 x 5e-4 */
void fig6() {
    const double W_A = 1.17e-3, G_A = -4.8e-9, G_B = -1.1e-8;
    const double KD = 6.3, DD = 7.6e-6, G_C = 2.8e-6;
    const double W_D1 = 1.09e-3, G_D1 = 1.4e-5 * 5e-4, G_D2 = 1.5e-5 * 5e-4;
    const double KP[2] = {0.85, 3.85}, DP[2] = {9.8e-8, 1.0e-4};
    const double T_END = 4.5e9;

    const double w_d2 = W_A + DP[0] - W_D1;
    const double w_b  = W_D1 + w_d2 + DP[1];
    const double w_c  = W_A + w_b - DD;

    // Initial energies read off the paper: a = b = d1 = 1e-9, c = d2 = 1e-11.
    const double q0 = std::sqrt(1e-9), q1 = std::sqrt(1e-11);

    // (a) direct triplet only -- unstable, as in Figure 1(b)
    panel("fig6a", Network({{"a", -W_A, G_A}, {"b", -w_b, G_B}, {"c", w_c, G_C}},
                           {Triplet{2, 0, 1, KD}}),
          State{q0, q0, q1}, T_END, 20000, 1e-10, 1e0, {}, DAY);

    // (b) five-mode mixed network
    const Network net({{"a", -W_A, G_A}, {"b", -w_b, G_B}, {"c", w_c, G_C},
                       {"d1", W_D1, G_D1}, {"d2", w_d2, G_D2}},
                      {Triplet{2, 0, 1, KD},
                       Triplet{0, 3, 4, KP[0]},
                       Triplet{1, 3, 4, KP[1]}});
    for (int i = 0; i < 2; ++i)
        line("6b", i ? "E_th_b" : "E_th_a",
             stab::threshold_energy(KP[i], W_D1, w_d2, G_D1, G_D2, DP[i]));
    panel("fig6b", net, State{q0, q0, q1, q0, q1}, T_END, 5000, 1e-10, 1e0, {}, DAY);
    std::printf("  mu = %.0f\n", stab::mu(KD, w_c, DD, G_C));
}

}  // namespace

int main(int argc, char** argv) {
    CLI::App app{"MW25 Figures 1-6, integrated with the C++ network"};
    std::string out = "out/mw25_cpp";
    std::vector<std::string> only;
    app.add_option("--out", out);
    app.add_option("--only", only, "figure numbers, e.g. --only 3 6")->expected(-1);
    CLI11_PARSE(app, argc, argv);

    OUT = out;
    std::filesystem::create_directories(OUT);
    auto want = [&](const char* n) {
        return only.empty() || std::find(only.begin(), only.end(), n) != only.end();
    };

    std::printf("MW25 figures -> %s\n", OUT.c_str());
    if (want("1")) fig1();
    if (want("2")) fig2();
    if (want("3")) fig3();
    if (want("4")) fig4();
    if (want("5")) fig5();
    if (want("6")) fig6();

    /* Merge rather than truncate: a --only run regenerates some panels and
       must not drop the lines belonging to the others. */
    const std::filesystem::path lp = OUT / "lines.csv";
    if (std::filesystem::exists(lp)) {
        const csv::Table old = csv::read(lp);
        for (size_t i = 0; i < old.rows(); ++i) {
            const std::string& pan = old.text("panel")[i];
            const std::string& lab = old.text("label")[i];
            bool superseded = false;
            for (const auto& l : LINES) superseded |= l[0] == pan && l[1] == lab;
            if (!superseded)
                LINES.push_back({pan, lab, old.text("value")[i]});
        }
    }
    std::sort(LINES.begin(), LINES.end());
    csv::Writer w(lp, {"panel", "label", "value"});
    for (const auto& l : LINES) w.row(l[0], l[1], l[2]);
    std::printf("  lines.csv  %zu entries\n", LINES.size());
    return 0;
}
