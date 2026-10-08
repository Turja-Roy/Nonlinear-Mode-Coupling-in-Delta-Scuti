/* Integrate an arbitrary network read from a plain-text file:

       ./network_run data/fig6_network.data --out out/network_run
       python3 ../scripts/network_plot.py --data out/network_run

   File format (see network_build.cpp's write_network): lines starting with
   '#' are comments, ignored, and so is the first data line (a norm-type
   label, not numbers). Then the mode table, one mode per line:

       id  gen  n  l  m  omega  gamma  [flin]

   (ids need not be 0-based or contiguous; gen and flin are not used here),
   then the triplet table:

       id1  id2  id3  kappa  [Delta/omega_1]

   A line's field count picks the table -- 7 or 8 is a mode, 4 or 5 a
   triplet. id2 == id3 is a literal self-coupling (id1 -> id2 + id2).

   omega sign: taken from the file. network_build writes each mode twice,
   +omega and -omega (q- = conj(q+), amplitude.hpp), paired here by (n, l, m);
   the - copy starts at the conjugate of the + one and outputs list the +
   copy only, under its id; energies_all.csv has every copy. A file with no negative omega at all (the old
   data files) is signed the MW25 way instead, negative iff gamma is. */

#include "amplitude.hpp"
#include "csv.hpp"

#include <CLI/CLI.hpp>

#include <algorithm>
#include <array>
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
    std::vector<int> mirror;                     // index of the -omega copy, -1 if none
    std::vector<int> shown;                      // the modes the outputs list
};

NetworkFile read_network (const std::filesystem::path& p) {
    std::ifstream f(p);
    if (!f) throw std::runtime_error("network_run: cannot open " + p.string());

    NetworkFile net;
    std::map<long, int> index;                   // file id -> Network mode index
    std::vector<std::array<double, 3>> nlm;
    bool skipped_norm_type = false;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (!skipped_norm_type) { skipped_norm_type = true; continue; }

        std::istringstream ss(line);
        std::vector<double> field;
        for (double v; ss >> v; ) field.push_back(v);

        if (field.size() == 7 || field.size() == 8) {         // id gen n l m omega gamma [flin]
            index[long(field[0])] = int(net.modes.size());
            nlm.push_back({field[2], field[3], field[4]});
            net.modes.push_back(amp::Mode{std::to_string(long(field[0])), field[5], field[6]});
        } else if (field.size() == 4 || field.size() == 5) {  // id1 id2 id3 kappa [Delta/w]
            net.triplets.push_back(Triplet{index.at(long(field[0])), index.at(long(field[1])),
                                           index.at(long(field[2])), field[3]});
        }
    }

    const bool doubled = std::any_of(net.modes.begin(), net.modes.end(),
                                     [](const amp::Mode& m) { return m.omega < 0.0; });
    if (!doubled)
        for (auto& m : net.modes) if (m.gamma < 0.0) m.omega = -m.omega;
    net.mirror.assign(net.modes.size(), -1);
    std::map<std::array<double, 3>, int> minus;
    if (doubled)
        for (size_t i = 0; i < net.modes.size(); ++i)
            if (net.modes[i].omega < 0.0) minus[nlm[i]] = int(i);
    for (size_t i = 0; i < net.modes.size(); ++i) {
        if (doubled && net.modes[i].omega < 0.0) {
            if (!minus.count(nlm[i]) || minus.at(nlm[i]) != int(i))
                throw std::runtime_error("network_run: two -omega copies of one (n, l, m)");
            continue;
        }
        net.shown.push_back(int(i));
        if (doubled) {
            const auto it = minus.find(nlm[i]);
            if (it == minus.end()) throw std::runtime_error("network_run: unpaired +omega mode");
            net.mirror[i] = it->second;
        }
    }
    if (doubled && net.shown.size() * 2 != net.modes.size())
        throw std::runtime_error("network_run: unpaired -omega mode");
    return net;
}

}  // namespace

int main (int argc, char** argv) {
    CLI::App app{"Integrate an arbitrary network read from a plain-text file"};
    std::string src, out = "out/network_run";
    double q0 = 1e-6, q0_d = 0.0, t_end = 0.0, rtol = 1e-10, e_max = 0.0;
    int n_out = 4000;
    unsigned seed = 0;
    app.add_option("file", src, "network file, see data/fake_net_101.data")->required();
    app.add_option("--out", out);
    app.add_option("--q0", q0, "initial |q| of the self-excited modes");
    app.add_option("--q0-daughter", q0_d, "initial |q| of the damped modes; default --q0");
    app.add_option("--seed", seed, "nonzero: random initial phases from this seed");
    app.add_option("--e-max", e_max, "stop once sum E crosses this [E_star]; 0 = never");
    app.add_option("--t-end", t_end, "integration length; default 20 / min|gamma|");
    app.add_option("--n-out", n_out);
    app.add_option("--rtol", rtol);
    CLI11_PARSE(app, argc, argv);

    const NetworkFile nf = read_network(src);
    const Network net(nf.modes, nf.triplets);
    const size_t copies = nf.modes.size() / nf.shown.size();       // 2 for a doubled file
    const long n_modes = long(nf.shown.size()), n_trip = long(nf.triplets.size() / copies);
    std::printf("%s: %ld modes, %ld triplets%s\n", src.c_str(), n_modes, n_trip,
                copies == 2 ? " (each as +-omega)" : "");

    if (t_end == 0.0) {
        double min_g = std::numeric_limits<double>::infinity();
        for (const auto& m : net.modes()) min_g = std::min(min_g, std::abs(m.gamma));
        t_end = 20.0 / min_g;
    }

    amp::Options opt;
    opt.n_out = n_out;
    opt.rtol = rtol;
    opt.e_max = e_max * double(copies);         // the - copies hold as much again
    if (q0_d == 0.0) q0_d = q0;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> phase(0.0, 2.0 * M_PI);
    State y0(size_t(net.size()));
    for (int i : nf.shown) {
        y0[size_t(i)] = std::polar(net.modes()[size_t(i)].gamma < 0.0 ? q0 : q0_d,
                                   seed ? phase(rng) : 0.0);
        if (nf.mirror[size_t(i)] >= 0) y0[size_t(nf.mirror[size_t(i)])] = std::conj(y0[size_t(i)]);
    }
    const amp::Solution sol = net.integrate(y0, t_end, opt);

    std::filesystem::create_directories(out);
    std::vector<std::string> hdr = {"t"};
    for (int i : nf.shown) hdr.push_back("E_" + net.modes()[size_t(i)].name);
    csv::Writer w(std::filesystem::path(out) / "energies.csv", hdr);
    for (size_t j=0 ; j<sol.t.size() ; j++) {
        const eig::ArrayXd E = net.energy(sol.y[j]);
        std::string row = csv::fmt(sol.t[j]);
        for (int i : nf.shown) row += "," + csv::fmt(E[i]);
        w.row(row);
    }
    std::printf("-> %s/energies.csv (t_end %.3e, %zu rows)\n", out.c_str(), t_end, sol.t.size());

    if (copies == 2) {
        hdr = {"t"};
        for (const auto& m : net.modes()) hdr.push_back("E_" + m.name);
        csv::Writer wa(std::filesystem::path(out) / "energies_all.csv", hdr);
        for (size_t j=0 ; j<sol.t.size() ; j++) {
            const eig::ArrayXd E = net.energy(sol.y[j]);
            std::string row = csv::fmt(sol.t[j]);
            for (eig::Index i=0 ; i<E.size() ; i++) row += "," + csv::fmt(E[i]);
            wa.row(row);
        }
    }

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
        for (int i : nf.shown) if (net.modes()[size_t(i)].gamma < 0.0) p += E[i];
        E_par.push_back(p);
    }
    double P = 0.0, D = 0.0, par = 0.0, top_d = 0.0;
    for (int i : nf.shown) {
        const double g = net.modes()[size_t(i)].gamma;
        (g < 0.0 ? P : D) += 2.0 * std::abs(g) * mean[i];
        if (g < 0.0) par += mean[i];
        else top_d = std::max(top_d, mean[i]);
    }
    long active = 0;
    for (int i : nf.shown)
        active += net.modes()[size_t(i)].gamma >= 0.0 && mean[i] > 1e-3 * top_d;
    const auto mm = std::minmax_element(E_par.begin(), E_par.end());
    // Stopped early: energy cap crossed, or the amplitudes overflowed.
    const bool runaway = sol.t.back() < t_end * (1.0 - 0.5 / (n_out - 1));
    const std::string stop = runaway ? "runaway" : "end";
    const double swing = (*mm.second - *mm.first) / par;
    std::printf("  <E_parents> %.3e E_star, swing %.3f, D/P %.3f, %ld active daughters, "
                "stop: %s at t %.3e\n", par, swing, D / P, active, stop.c_str(), sol.t.back());
    csv::Writer s(std::filesystem::path(out) / "summary.csv",
                  {"file", "n_modes", "n_triplets", "t_end", "q0", "q0_daughter", "seed",
                   "runaway", "stop", "E_parents", "swing", "driving", "dissipation", "n_active"});
    s.row(src, n_modes, n_trip, sol.t.back(), q0, q0_d,
          long(seed), runaway, stop, par, swing, P, D, active);
    return 0;
}
