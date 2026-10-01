/* Integrate an arbitrary network read from a plain-text file:

       ./network_run data/fig6_network.data --out out/network_run
       python3 ../scripts/network_plot.py --data out/network_run

   File format (see data/fake_net_101.data): lines starting with '#' are
   comments, ignored, and so is the first data line (a norm-type label, not
   numbers). Then the mode table, one mode per line:

       id  gen  n  l  m  omega  gamma  flin

   (ids need not be 0-based or contiguous; gen and flin are not used here),
   then the triplet table:

       id1  id2  id3  kappa

   A line's field count picks the table -- 8 is a mode, 4 is a triplet -- so
   no comment-block bookkeeping is needed beyond skipping '#' lines. id2 ==
   id3 in a triplet row is a literal self-coupling (id1 -> id2 + id2) and
   amp::Network handles it as such, not as an error.

   omega/gamma sign: the file's omega is a plain magnitude; the amplitude
   equations need it signed by the mode's own stability role, negative iff
   gamma is (self-excited), exactly as amp::from_triplets does it for the
   stellar-model pipeline -- see amplitude.hpp's header comment for why. */

#include "amplitude.hpp"
#include "csv.hpp"

#include <CLI/CLI.hpp>

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

        if (field.size() == 8) {                  // id gen n l m omega gamma flin
            const long id = long(field[0]);
            const double gamma = field[6];
            const double omega = (gamma < 0.0 ? -1.0 : 1.0) * std::abs(field[5]);
            index[id] = int(net.modes.size());
            net.modes.push_back(amp::Mode{std::to_string(id), omega, gamma});
        } else if (field.size() == 4) {            // id1 id2 id3 kappa
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
    std::printf("%s: %zu modes, %zu triplets\n", src.c_str(), nf.modes.size(), nf.triplets.size());

    if (t_end == 0.0) {
        double min_g = std::numeric_limits<double>::infinity();
        for (const auto& m : net.modes()) min_g = std::min(min_g, std::abs(m.gamma));
        t_end = 20.0 / min_g;
    }

    amp::Options opt;
    opt.n_out = n_out;
    opt.rtol = rtol;
    opt.e_max = e_max;
    if (q0_d == 0.0) q0_d = q0;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> phase(0.0, 2.0 * M_PI);
    State y0(size_t(net.size()));
    for (size_t i = 0; i < y0.size(); ++i)
        y0[i] = std::polar(net.modes()[i].gamma < 0.0 ? q0 : q0_d, seed ? phase(rng) : 0.0);
    const amp::Solution sol = net.integrate(y0, t_end, opt);

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
        for (int i = 0; i < net.size(); ++i) if (net.modes()[i].gamma < 0.0) p += E[i];
        E_par.push_back(p);
    }
    double P = 0.0, D = 0.0, par = 0.0, top_d = 0.0;
    for (int i = 0; i < net.size(); ++i) {
        const double g = net.modes()[i].gamma;
        (g < 0.0 ? P : D) += 2.0 * std::abs(g) * mean[i];
        if (g < 0.0) par += mean[i];
        else top_d = std::max(top_d, mean[i]);
    }
    long active = 0;
    for (int i = 0; i < net.size(); ++i)
        active += net.modes()[i].gamma >= 0.0 && mean[i] > 1e-3 * top_d;
    const auto mm = std::minmax_element(E_par.begin(), E_par.end());
    const bool runaway = sol.t.back() < t_end;
    const double swing = (*mm.second - *mm.first) / par;
    std::printf("  <E_parents> %.3e E_star, swing %.3f, D/P %.3f, %ld active daughters%s\n",
                par, swing, D / P, active, runaway ? ", RUNAWAY" : "");
    csv::Writer s(std::filesystem::path(out) / "summary.csv",
                  {"file", "n_modes", "n_triplets", "t_end", "q0", "q0_daughter", "seed",
                   "runaway", "E_parents", "swing", "driving", "dissipation", "n_active"});
    s.row(src, long(net.size()), long(net.triplets().size()), sol.t.back(), q0, q0_d,
          long(seed), runaway, par, swing, P, D, active);
    return 0;
}
