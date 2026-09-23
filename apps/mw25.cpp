/* Mourabit & Weinberg (2025), ApJ 986:33, Figures 1-6, integrated in C++.

       ./mw25 --out out/mw25_cpp
       python3 scripts/mw25_plot.py --data out/mw25_cpp

   This writes one CSV per panel (t and the mode energies) plus a single
   lines.csv carrying the horizontal thresholds and equilibria. Plotting stays
   in Python; the point of running it here is that the trajectories then come
   out of the same Network the pipeline uses.

   Sign conventions. MW25 write signed frequencies whose sum is the detuning,
   so a parent carries omega < 0. amp::Network stores |omega| with the sum mode
   first in every triplet, so Delta = w_b + w_c - w_a up to an overall sign
   that does not affect |q|. */

#include "amplitude.hpp"
#include "csv.hpp"
#include "stability.hpp"

#include <CLI/CLI.hpp>

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
        /* Seems like the MW25 fig1 is actually using w_a = w_b = 1.0, w_c = 2.0,
            not -0.5, -0.5, 1.0 as stated in the caption
            and \gamma_a = gamma_b = -0.01, gamma_c = 0.1, as opposed to
            -0.001, -0.001, 0.01 in the caption.
            However, unsure about the negative frequency there;
            paper mentions \omega_a = \omega_b = -0.5 */
        // Network net({ {"c", 1.0, 0.01}, {"a", 0.5, -0.001} }, {Triplet{0, 1, 1, 1.0}});
        Network net({ {"c", 2.0, 0.10}, {"a", 1.0, -0.010} }, {Triplet{0, 1, 1, 1.0}});
        const double d = net.detuning(net.triplets()[0]);
        // MW25 Eq. (6) for the sum-slot mode; its two partners are the parent
        // twice, and both are driven, which threshold_energy() refuses.
        const double Ec = (-0.001) * (-0.001) / (4 * 1.0 * 0.5 * 0.5)
                        * (1 + d * d / (0.008 * 0.008));
        const double Ea = (-0.001) * (0.01) / (4 * 1.0 * (-0.5) * 1.0)
                        * (1 + d * d / (0.008 * 0.008));
        line("1a", "E_c_eq", Ec);
        line("1a", "E_a_eq", Ea);
        panel("fig1a", net, State(2, 1e-3), 6000.0, 4000, 1e-11);
    }
    // (b) two distinct parents: the extra degree of freedom breaks it
    // Network net({{"c", 1.0, 0.01}, {"a", 0.500, -0.0011}, {"b", 0.502, -0.001}},
    //             {Triplet{0, 1, 2, 1.0}});
    /* Same here */
    Network net({{"c", 2.0, 0.10}, {"a", 1.0, -0.010}, {"b", 1.004, -0.010}},
                {Triplet{0, 1, 2, 1.0}});
    panel("fig1b", net, State(3, 1e-3), 6000.0, 4000, 1e-11, 1e2);
}

void fig2() {
    /* Deterministic 1% spread/jitter. Otherwise all the daughter are identical and the
        roundoff in the integrator picks one arbitrarily, lets it grow first, and then
        the energy sloshes back and forth among the N identical channels
        that gives a bouncy plot. */
    unsigned long s = 12345;
    auto jitter = [&s] {
        s = (1103515245UL * s + 12345UL) & 0x7fffffffUL;
        return 0.01 * (2.0 * (double(s) / 2147483647.0) - 1.0);
    };
    for (int N : {2, 10, 50}) {
        std::vector<amp::Mode> modes = {{"a", 0.5, -0.001}, {"b", 0.5, -0.0011}};
        std::vector<Triplet> trs;
        for (int i=0 ; i<N ; i++) {
            const double wd = 1.0 * (1 + jitter()), gd = 0.01 * (1 + jitter());
            const double kd = 1.0 * (1 + jitter());
            // const double wd = 1.0, gd = 0.01;
            // const double kd = 1.0;
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
constexpr double OMEGA[4] = {1.000, 1.001, 2.0, 0.5};
constexpr double GAMMA[4] = {-0.01, -0.01, 0.1, 0.1};

Network mixed(double k_direct, double k_param, const double g[4]) {
    return Network({{"a", OMEGA[0], g[0]}, {"b", OMEGA[1], g[1]},
                    {"c", OMEGA[2], g[2]}, {"d", OMEGA[3], g[3]}},
                   {Triplet{2, 0, 1, k_direct},      // a + b -> c
                    Triplet{0, 3, 3, k_param},       // a -> d + d
                    Triplet{1, 3, 3, k_param}});     // b -> d + d
}

void mixed_lines(const std::string& tag, const Network& net, double k_param,
                 const double g[4]) {
    const double dp = net.detuning(net.triplets()[1]);
    line(tag, "E_th", stab::threshold_energy(k_param, OMEGA[3], OMEGA[3], g[3], g[3], dp));
    line(tag, "E_d_eq",
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

/* --------------------------------------------------------------- Figure 6 */

/* MW25 Section 6, their 2.0 Msun, Teff = 7202 K, log g = 3.80 model.
   Parents a, b and direct daughter c: l = {3, 1, 2}, n = {0, 3, 8}.
   Parametric daughter pair: l = {14, 15}, n = {-282, -12}.

   The published frequencies are quoted to three significant figures and do not
   reproduce the published detunings: 1.17 + 1.27 - 2.42 = 2e-5, not the stated
   7.6e-5 rad/s. Delta drives the dynamics at leading order while the
   individual omega only enter as prefactors, so each triplet's sum-slot
   frequency is set here from the quoted Delta. The shifts are ~2%. */
// void fig6() {
//     constexpr double W_A = 1.17e-3, G_A = -4.8e-9, G_B = -1.1e-8, G_C = 2.8e-5;
//     constexpr double KD = 6.3, DD = 7.6e-5;
//     constexpr double W_D1 = 1.09e-3, G_D1 = 1.4e-5, G_D2 = 1.5e-5;
//     constexpr double KP[2] = {0.85, 3.85}, DP[2] = {9.8e-8, 1.0e-4};
// 
//     const double w_d2 = W_A + DP[0] - W_D1;         // 8.01e-5 vs the quoted 8.22e-5
//     const double w_b  = W_D1 + w_d2 + DP[1];
//     const double w_c  = W_A + w_b - DD;             // 2.364e-3 vs the quoted 2.42e-3
// 
//     // (a) direct triplet only -- unstable, as in Figure 1(b)
//     panel("fig6a", Network({{"a", W_A, G_A}, {"b", w_b, G_B}, {"c", w_c, G_C}},
//                            {Triplet{2, 0, 1, KD}}),
//           State(3, 1e-3), 3.0e9, 3000, 1e-8, 1e0, {}, DAY);
// 
//     // (b) five-mode mixed network. Seeded near the parametric threshold: the
//     // stretch below it is pure linear growth and costs more to integrate than
//     // the limit cycle it leads to.
//     const Network net({{"a", W_A, G_A}, {"b", w_b, G_B}, {"c", w_c, G_C},
//                        {"d1", W_D1, G_D1}, {"d2", w_d2, G_D2}},
//                       {Triplet{2, 0, 1, KD},
//                        Triplet{0, 3, 4, KP[0]},
//                        Triplet{1, 3, 4, KP[1]}});
//     for (int i = 0; i < 2; ++i)
//         line("6b", i ? "E_th_b" : "E_th_a",
//              stab::threshold_energy(KP[i], W_D1, w_d2, G_D1, G_D2, DP[i]));
//     panel("fig6b", net, State(5, 1e-3), 3.0e9, 4000, 1e-8, 1e0, {}, DAY);
// 
//     std::printf("\n  note: with the published gamma_a = %.1e 1/s the parent e-folds in"
//                 " %.0f d,\n        so these panels span ~%.0f d where MW25's axis runs"
//                 " to 2.8e8 d.\n        Their time axis and their quoted rates differ by"
//                 " ~1e4; the rates are used here.\n",
//                 G_A, 1 / std::abs(G_A) / DAY, 3.0e9 / DAY);
// }

/* Variants test which quoted number is the typo (published mu ~ 2000, but the
   quoted values give ~190):
     0  as published
     1  Delta_direct = 7.6e-6
     2  1 + gamma_c = 2.8e-6
     3  kappa_direct = 63 */
void fig6(int v, double t_end) {
    constexpr double W_A = 1.17e-3, G_A = -4.8e-9, G_B = -1.1e-8;
    constexpr double W_D1 = 1.09e-3, G_D1 = 1.4e-5, G_D2 = 1.5e-5;
    constexpr double KP[2] = {0.85, 3.85}, DP[2] = {9.8e-8, 1.0e-4};
    const double KD  = v == 3 ? 63.0 : 6.3;
    const double DD  = v == 1 || v == 2 ? 7.6e-6 : 7.6e-5;
    const double G_C = v == 2 ? 2.8e-6 : 2.8e-5;

    const double w_d2 = W_A + DP[0] - W_D1;
    const double w_b  = W_D1 + w_d2 + DP[1];
    const double w_c  = W_A + w_b - DD;

    // Seeded below the direct-coupling saddle (E ~ 5e-9) so the parent with
    // the larger omega and driving, b, is the one that runs away (MW25 fn. 1).
    panel("fig6a", Network({{"a", W_A, G_A}, {"b", w_b, G_B}, {"c", w_c, G_C}},
                           {Triplet{2, 0, 1, KD}}),
          State{1e-4, 1e-4, 3e-6}, t_end, 3000, 1e-8, 1e0, {}, DAY);

    const Network net({{"a", W_A, G_A}, {"b", w_b, G_B}, {"c", w_c, G_C},
                       {"d1", W_D1, G_D1}, {"d2", w_d2, G_D2}},
                      {Triplet{2, 0, 1, KD},
                       Triplet{0, 3, 4, KP[0]},
                       Triplet{1, 3, 4, KP[1]}});
    for (int i = 0; i < 2; ++i)
        line("6b", i ? "E_th_b" : "E_th_a",
             stab::threshold_energy(KP[i], W_D1, w_d2, G_D1, G_D2, DP[i]));
    panel("fig6b", net, State{1e-4, 1e-4, 3e-6, 1e-5, 1e-5}, t_end, 4000, 1e-8, 1e0, {}, DAY);
    std::printf("  variant %d: mu = %.0f\n", v, stab::mu(KD, w_c, DD, G_C));
}

}  // namespace

int main(int argc, char** argv) {
    CLI::App app{"MW25 Figures 1-6, integrated with the C++ network"};
    std::string out = "out/mw25_cpp";
    std::vector<std::string> only;
    app.add_option("--out", out);
    app.add_option("--only", only, "figure numbers, e.g. --only 3 6")->expected(-1);
    int v6 = 0;
    double t6 = 3.0e9;
    app.add_option("--v6", v6, "fig6 parameter variant 0-3");
    app.add_option("--t6", t6, "fig6 t_end [s]");
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
    if (want("6")) fig6(v6, t6);

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
