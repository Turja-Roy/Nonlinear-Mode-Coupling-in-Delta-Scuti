/* dq_dt against dA_dt: each case integrated as q (rtol 1e-13, 1e-12) and as A
   (rtol 1e-10, 1e-8), compared as complex q = A exp(-i w t).

       ./slow_amplitude --out out/slow_amplitude
       ./slow_amplitude --net data/atol_bench/1p2d_M2.0T_eth.data --out out/slow_amplitude
       python3 scripts/slow_amplitude_plot.py --data out/slow_amplitude
*/

#include "amplitude.hpp"
#include "csv.hpp"
#include "stability.hpp"

#include <CLI/CLI.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>

using amp::Network;
using amp::State;
using amp::Triplet;

namespace {

struct Form { const char* name; bool slow; double rtol; };
constexpr Form FORMS[4] = {{"q13", false, 1e-13}, {"q12", false, 1e-12}, {"A10", true, 1e-10}, {"A8", true, 1e-8}};
constexpr int PAIRS[4][2] = {{0, 2}, {0, 3}, {1, 2}, {1, 0}};

struct Case {
    std::string name;
    Network net;
    std::vector<double> amp0;
    double t_end, e_max = 0.0;
};

// parents at q0, daughters 1e-3 below in energy
std::vector<double> parent_daughter (const Network& net, double q0) {
    std::vector<double> a;
    for (const auto& m : net.modes()) a.push_back(m.gamma < 0.0 ? q0 : std::sqrt(1e-3) * q0);
    return a;
}

std::vector<Case> toys () {
    std::vector<Case> cs;
    auto three = [](double delta, const char* name) {
        const Network n = amp::three_mode({-2.0 + delta, 1.1, 0.9}, {-0.001, 0.01, 0.01}, 1.0);
        return Case{name, n, parent_daughter(n, 1e-3), 2e4};
    };
    // MW25 sec. 5.1: steady for |Delta| > gamma_d, else a cycle
    cs.push_back(three(0.05, "T1_parametric_steady"));
    cs.push_back(three(0.002, "T2_parametric_cycle"));
    const Network sc = amp::self_coupled(-2.0 + 0.05, 1.0, -0.001, 0.01, 1.0);
    cs.push_back({"T3_self_coupled", sc, parent_daughter(sc, 1e-3), 2e4});
    cs.push_back({"T4_mw25_fig1a", Network({{"c", 2.0, 0.1}, {"a", -1.0, -0.01}}, {Triplet{0, 1, 1, 1.0}}),
                  {1e-3, 1e-3}, 6000.0});
    cs.push_back({"T5_mw25_fig1b",
                  Network({{"c", 2.0, 0.1}, {"a", -1.0, -0.01}, {"b", -1.004, -0.01}}, {Triplet{0, 1, 2, 1.0}}),
                  {1e-3, 1e-3, 1e-3}, 6000.0, 1e2});
    cs.push_back({"T6_mw25_fig3",
                  Network({{"a", -1.000, -0.01}, {"b", -1.001, -0.01}, {"c", 2.0, 0.1}, {"d", 0.5, 0.1}},
                          {Triplet{2, 0, 1, 1.0}, Triplet{0, 3, 3, 1.0}, Triplet{1, 3, 3, 1.0}}),
                  {1e-3, 1e-3, 1e-3, 1e-3}, 2000.0});
    return cs;
}

void run (const Case& c, const std::filesystem::path& out, int n_out, double max_s) {
    const Network& net = c.net;
    const int n = net.size();
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> phase(0.0, 2.0 * M_PI);
    State y0(n);
    for (int i=0 ; i<n ; i++) y0[i] = std::polar(c.amp0[i], phase(rng));

    std::vector<amp::Solution> sol;
    std::vector<double> secs;
    for (const Form& f : FORMS) {
        amp::Options opt;
        opt.n_out = n_out;
        opt.rtol = f.rtol;
        opt.slow = f.slow;
        opt.e_max = c.e_max;
        opt.max_seconds = max_s;
        const auto t0 = std::chrono::steady_clock::now();
        sol.push_back(net.integrate(y0, c.t_end, opt));
        secs.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        std::printf("  %-22s %-3s: %6zu outputs to t %.4g, %ld steps, %.1f s\n", c.name.c_str(), f.name,
                    sol.back().t.size(), sol.back().t.back(), sol.back().steps, secs.back());
    }

    // e_max puts each run's last output at its own crossing
    size_t n_rows = sol[0].t.size();
    for (const auto& s : sol) n_rows = std::min(n_rows, s.t.size());
    std::vector<size_t> rows;
    for (size_t j=0 ; j<n_rows ; j++)
        if (std::all_of(sol.begin(), sol.end(), [&](const amp::Solution& s) { return s.t[j] == sol[0].t[j]; }))
            rows.push_back(j);

    std::vector<std::string> hdr = {"t"};
    for (const Form& f : FORMS)
        for (const auto& m : net.modes()) hdr.push_back(std::string("E_") + f.name + "_" + m.name);
    for (const auto& p : PAIRS)
        for (const auto& m : net.modes())
            hdr.push_back(std::string("d_") + FORMS[p[0]].name + "_" + FORMS[p[1]].name + "_" + m.name);
    csv::Writer w(out / (c.name + ".csv"), hdr);
    std::vector<double> d_max(std::size(PAIRS), 0.0);
    for (size_t j : rows) {
        const double t = sol[0].t[j];
        std::string row = csv::fmt(t);
        std::vector<State> q;
        for (int k=0 ; k<4 ; k++) {
            q.push_back(FORMS[k].slow ? net.to_q(t, sol[k].y[j]) : sol[k].y[j]);
            for (int i=0 ; i<n ; i++) row += "," + csv::fmt(std::norm(q.back()[i]));
        }
        for (size_t p=0 ; p<std::size(PAIRS) ; p++)
            for (int i=0 ; i<n ; i++) {
                const State &x = q[PAIRS[p][0]], &y = q[PAIRS[p][1]];
                const double d = std::abs(x[i] - y[i]) / std::abs(x[i]);
                d_max[p] = std::max(d_max[p], d);
                row += "," + csv::fmt(d);
            }
        w.row(row);
    }

    std::vector<double> E_fix(size_t(n), std::nan(""));
    if (net.triplets().size() == 1) {
        const Triplet& tr = net.triplets()[0];
        const auto& md = net.modes();
        const auto E = stab::equilibrium_energies(tr.kappa, {md[tr.a].omega, md[tr.b].omega, md[tr.c].omega},
                                                  {md[tr.a].gamma, md[tr.b].gamma, md[tr.c].gamma},
                                                  net.detuning(tr), 1.0, tr.b == tr.c);
        E_fix[tr.a] = E[0];
        E_fix[tr.b] = E[1];
        E_fix[tr.c] = E[2];
    }
    csv::Writer wm(out / (c.name + "_modes.csv"), {"mode", "omega", "gamma", "E_fixed_point"});
    for (int i=0 ; i<n ; i++) wm.row(net.modes()[i].name, net.modes()[i].omega, net.modes()[i].gamma, E_fix[i]);

    csv::Writer wr(out / (c.name + "_runs.csv"), {"run", "slow", "rtol", "t_end", "t_last", "stop", "rows_compared",
                                                 "seconds", "steps", "rhs_calls", "failed_steps"});
    for (int k=0 ; k<4 ; k++) {
        const amp::Solution& s = sol[k];
        const std::string stop = s.failed ? "failed" : s.out_of_time ? "wall"
                                 : s.t.back() < c.t_end * (1.0 - 0.5 / (n_out - 1)) ? "e_max" : "end";
        wr.row(std::string(FORMS[k].name), FORMS[k].slow, FORMS[k].rtol, c.t_end, s.t.back(), stop,
               long(rows.size()), secs[k], s.steps, s.rhs_calls, s.failed_steps);
    }
    csv::Writer wp(out / (c.name + "_pairs.csv"), {"x", "y", "d_max"});
    for (size_t p=0 ; p<std::size(PAIRS) ; p++) {
        const char *x = FORMS[PAIRS[p][0]].name, *y = FORMS[PAIRS[p][1]].name;
        wp.row(std::string(x), std::string(y), d_max[p]);
        std::printf("  %-22s max |q_%s - q_%s| / |q_%s|: %.2e\n", c.name.c_str(), x, y, x, d_max[p]);
    }
}

}  // namespace

int main (int argc, char** argv) {
    CLI::App app{"dq_dt against dA_dt"};
    std::string out = "out/slow_amplitude";
    std::vector<std::string> nets;
    double t_end = 1e12, max_s = 10800.0;
    int n_out = 4000;
    app.add_option("--out", out, "default out/slow_amplitude");
    app.add_option("--net", nets, "network files instead of the toys");
    app.add_option("--t-end", t_end, "network files [s], default 1e12");
    app.add_option("--n-out", n_out, "default 4000");
    app.add_option("--max-seconds", max_s, "wall-time cap per run, default 10800; 0 = none");
    CLI11_PARSE(app, argc, argv);

    std::filesystem::create_directories(out);
    std::vector<Case> cs;
    if (nets.empty()) cs = toys();
    for (const auto& f : nets) {
        const Network net = amp::read_network(f);
        cs.push_back({std::filesystem::path(f).stem().string(), net, parent_daughter(net, 1e-6), t_end});
    }
    for (const Case& c : cs) run(c, out, n_out, max_s);
    return 0;
}
