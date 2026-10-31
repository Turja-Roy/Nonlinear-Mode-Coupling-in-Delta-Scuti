/* Integrate an arbitrary network read from a plain-text file:

       ./network_run data/fig6_network.data --out out/network_run
       python3 ../scripts/network_plot.py --data out/network_run

   File format (see network_build.cpp's write_network): lines starting with
   '#' are comments, ignored, and so is the first data line (a norm-type
   label, not numbers). Then the mode table, one mode per line:

       id  gen  n  l  m  omega  gamma  [flin]

   (ids need not be 0-based or contiguous; flin is not used here),
   then the triplet table:

       id1  id2  id3  kappa  [Delta/omega_1]

   A line's field count picks the table -- 7 or 8 is a mode, 4 or 5 a
   triplet. id2 == id3 is a literal self-coupling (id1 -> id2 + id2).

   omega sign: by generation, -|omega| for even gen, +|omega| for odd
   (amplitude.hpp); the file's omega is a magnitude. A negative one means a
   +-omega file from an older network_build, which is refused.

   Starting energies: every parent at |q| = --q0, every daughter r times
   below in energy, E_d(0) = r E_p(0), r = --daughter-ratio <= 1, so no
   daughter starts above a parent. Each parent's E_th (stab::threshold_energy,
   lowest over its parametric pairs in the file) is logged beside its E(0). */

#include "amplitude.hpp"
#include "csv.hpp"
#include "stability.hpp"

#include <CLI/CLI.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <random>
#include <sstream>

using amp::Network;
using amp::State;
using amp::Triplet;

namespace {

struct NetworkFile {
    std::vector<amp::Mode> modes;
    std::vector<Triplet> triplets;
};

NetworkFile read_network (const std::filesystem::path& p) {
    std::ifstream f(p);
    if (!f) throw std::runtime_error("network_run: cannot open " + p.string());

    NetworkFile net;
    std::map<long, int> index;                   // file id -> Network mode index
    bool skipped_norm_type = false;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (!skipped_norm_type) { skipped_norm_type = true; continue; }

        std::istringstream ss(line);
        std::vector<double> field;
        for (double v; ss >> v; ) field.push_back(v);

        if (field.size() == 7 || field.size() == 8) {         // id gen n l m omega gamma [flin]
            if (field[5] < 0.0)
                throw std::runtime_error("network_run: negative omega in " + p.string()
                                         + ", a +-omega file from an older network_build; rebuild it");
            index[long(field[0])] = int(net.modes.size());
            const double sign = long(field[1]) % 2 == 0 ? -1.0 : 1.0;
            net.modes.push_back(amp::Mode{std::to_string(long(field[0])), sign * field[5], field[6]});
        } else if (field.size() == 4 || field.size() == 5) {  // id1 id2 id3 kappa [Delta/w]
            net.triplets.push_back(Triplet{index.at(long(field[0])), index.at(long(field[1])),
                                           index.at(long(field[2])), field[3]});
        }
    }

    return net;
}

}  // namespace

int main (int argc, char** argv) {
    CLI::App app{"Integrate an arbitrary network read from a plain-text file"};
    std::string src, out = "out/network_run";
    double q0 = 1e-6, ratio = 1e-3, t_end = 0.0, rtol = 1e-10, atol = 1e-300,
           e_max = 0.0, max_s = 0.0;
    int n_out = 4000;
    unsigned seed = 0;
    app.add_option("file", src, "network file, see data/fake_net_101.data")->required();
    app.add_option("--out", out);
    app.add_option("--q0", q0, "initial |q| of every parent (self-excited mode)");
    app.add_option("--daughter-ratio", ratio, "E_daughter(0) / E_parent(0), at most 1")
        ->check(CLI::Range(0.0, 1.0));
    app.add_option("--seed", seed, "nonzero: random initial phases from this seed");
    app.add_option("--e-max", e_max, "stop once sum E crosses this [E_star]; 0 = never");
    app.add_option("--t-end", t_end, "integration length; default 20 / min|gamma|");
    app.add_option("--n-out", n_out);
    app.add_option("--rtol", rtol);
    app.add_option("--atol", atol, "absolute tolerance on |A| (amplitude.hpp)");
    app.add_option("--max-seconds", max_s, "stop integrating after this wall time, keeping outputs; 0 = never");
    CLI11_PARSE(app, argc, argv);

    const NetworkFile nf = read_network(src);
    const Network net(nf.modes, nf.triplets);
    const long n_modes = long(net.size()), n_trip = long(net.triplets().size());
    std::printf("%s: %ld modes, %ld triplets\n", src.c_str(), n_modes, n_trip);

    if (t_end == 0.0) {
        double min_g = std::numeric_limits<double>::infinity();
        for (const auto& m : net.modes()) min_g = std::min(min_g, std::abs(m.gamma));
        t_end = 20.0 / min_g;
    }

    amp::Options opt;
    opt.n_out = n_out;
    opt.rtol = rtol;
    opt.atol = atol;
    opt.max_seconds = max_s;
    opt.e_max = e_max;
    const auto& md = net.modes();
    const double inf = std::numeric_limits<double>::infinity();
    const double q0_d = std::sqrt(ratio) * q0;
    std::vector<double> e_th(md.size(), inf);
    for (const Triplet& tr : net.triplets()) {
        const double e = stab::threshold_energy(tr.kappa, md[tr.b].omega, md[tr.c].omega, md[tr.b].gamma,
                                                md[tr.c].gamma, net.detuning(tr));
        if (md[tr.a].gamma < 0.0 && e > 0.0) e_th[tr.a] = std::min(e_th[tr.a], e);
    }
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> phase(0.0, 2.0 * M_PI);
    State y0(size_t(net.size()));
    // phases on the +|omega| representation: a -omega mode starts at the conjugate
    for (size_t i = 0; i < y0.size(); ++i) {
        const bool parent = md[i].gamma < 0.0;
        const double q = parent ? q0 : q0_d;
        if (parent) std::printf("  parent %s: E_th %.3e, E(0) %.3e\n", md[i].name.c_str(), e_th[i], q * q);
        y0[i] = std::polar(q, (md[i].omega < 0.0 ? -1.0 : 1.0) * (seed ? phase(rng) : 0.0));
    }
    const auto clock0 = std::chrono::steady_clock::now();
    amp::Solution sol = net.integrate(y0, t_end, opt);
    // Stopped early by the energy cap (a solver failure or NaN is sol.failed instead).
    const bool runaway = !sol.out_of_time && !sol.failed && sol.t.back() < t_end * (1.0 - 0.5 / (n_out - 1));
    /* A fast parent can cross the cap within a few outputs, leaving a plot of
       two points; integrate again up to the crossing to resolve the growth.
       The counters below add up both passes: the cost of the run. */
    if (runaway && sol.t.size() < size_t(n_out / 4) && sol.t.back() > 0.0) {
        t_end = sol.t.back();
        if (max_s > 0.0)                         // one budget for both passes
            opt.max_seconds = std::max(1.0, max_s - std::chrono::duration<double>(
                                                std::chrono::steady_clock::now() - clock0).count());
        const amp::Solution first = sol;
        sol = net.integrate(y0, t_end, opt);
        sol.steps += first.steps; sol.rhs_calls += first.rhs_calls; sol.failed_steps += first.failed_steps;
        sol.nonlin_iters += first.nonlin_iters;
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count();
    std::printf("  integration %.1f s: %ld steps, %ld failed, %ld rhs calls, %ld nonlinear iterations\n",
                seconds, sol.steps, sol.failed_steps, sol.rhs_calls, sol.nonlin_iters);

    std::filesystem::create_directories(out);
    std::vector<std::string> hdr = {"t"};
    for (const auto& m : net.modes()) hdr.push_back("E_" + m.name);
    csv::Writer w(std::filesystem::path(out) / "energies.csv", hdr);
    for (size_t j=0 ; j<sol.t.size() ; j++) {
        const eig::ArrayXd E = net.energy(sol.y[j]);
        std::string row = csv::fmt(sol.t[j]);
        for (eig::Index i=0 ; i<E.size() ; i++) row += "," + csv::fmt(E[i]);
        w.row(row);
    }
    std::printf("-> %s/energies.csv (t_end %.3e, %zu rows)\n", out.c_str(), t_end, sol.t.size());

    /* Second-half averages. At a bounded state the parents' kappa-mechanism
       input equals what the damped modes dissipate, so D/P -> 1 is the check
       that the run has actually settled. */
    const size_t j0 = sol.t.size() / 2, nj = sol.t.size() - j0;
    eig::ArrayXd mean = eig::ArrayXd::Zero(net.size());
    std::vector<double> E_par;
    for (size_t j = j0; j < sol.t.size(); ++j) {
        const eig::ArrayXd E = net.energy(sol.y[j]);
        mean += E / double(nj);
        double p = 0.0;
        for (int i = 0; i < net.size(); ++i) if (net.modes()[size_t(i)].gamma < 0.0) p += E[i];
        E_par.push_back(p);
    }
    double P = 0.0, D = 0.0, par = 0.0, top_d = 0.0;
    for (int i = 0; i < net.size(); ++i) {
        const double g = net.modes()[size_t(i)].gamma;
        (g < 0.0 ? P : D) += 2.0 * std::abs(g) * mean[i];
        if (g < 0.0) par += mean[i];
        else top_d = std::max(top_d, mean[i]);
    }
    long active = 0;
    for (int i = 0; i < net.size(); ++i)
        active += net.modes()[size_t(i)].gamma >= 0.0 && mean[i] > 1e-3 * top_d;
    const auto mm = std::minmax_element(E_par.begin(), E_par.end());
    const std::string stop = sol.failed ? "failed" : sol.out_of_time ? "wall" : runaway ? "runaway" : "end";
    const double swing = (*mm.second - *mm.first) / par;
    std::printf("  <E_parents> %.3e E_star, swing %.3f, D/P %.3f, %ld active daughters, "
                "stop: %s at t %.3e\n", par, swing, D / P, active, stop.c_str(), sol.t.back());
    csv::Writer s(std::filesystem::path(out) / "summary.csv",
                  {"file", "n_modes", "n_triplets", "t_end", "q0", "q0_daughter", "daughter_ratio",
                   "seed", "runaway", "stop", "E_parents", "swing", "driving", "dissipation", "n_active",
                   "seconds", "steps", "failed_steps", "rhs_calls", "nonlin_iters"});
    s.row(src, n_modes, n_trip, sol.t.back(), q0, q0_d, ratio,
          long(seed), runaway, stop, par, swing, P, D, active,
          seconds, sol.steps, sol.failed_steps, sol.rhs_calls, sol.nonlin_iters);
    return 0;
}
