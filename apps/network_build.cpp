/* Parametric networks around one or two parents, written for network_run:

       ./network_build --parent 6,3 --parent 6,4 --rule eth delta -N 2 4 8 --out data/net
       ./network_run data/net/eth_cut0.15_N8.data --out out/net/eth_cut0.15_N8

   One file per (cut, rule, N) under --out, from one load of the model: the
   load is minutes, a network is milliseconds.

   Candidates are the triplets a parent heads, parent -> b + c with both
   daughters damped and |Delta| < cut. `rule` ranks them and pairs are taken in
   order until N daughters are in (one more if the last pair brings two):

       eth      E_th ascending -- kappa, gamma and Delta together
       delta    |Delta| ascending
       dgamma   |Delta| / (gamma_b + gamma_c) ascending, E_th without kappa
       random   uniform; the control for how much the ranking matters

   Closure then adds every triplet among the chosen modes with |Delta| below
   the closure cut, whichever leg selected them: the equations do not know why
   a mode is in the network. --no-closure keeps the selecting legs only. */

#include "csv.hpp"
#include "stability.hpp"

#include <CLI/CLI.hpp>

#include <cstdio>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>

using stab::Leg;

namespace {

Key parse_key (const std::string& s) {
    const size_t c = s.find(',');
    if (c == std::string::npos) throw std::runtime_error("--parent wants l,n: " + s);
    return {std::stoi(s.substr(0, c)), std::stoi(s.substr(c + 1))};
}

double kappa_of (const stab::KappaCache& kap, const RadialTriplet& t) {
    const auto k = t.keys();
    return kap.at(stab::leg_of(k[0], k[1], k[2])).first;
}

double e_th (const ModeMap& efs, const RadialTriplet& t, double kappa) {
    if (kappa == 0.0) return std::numeric_limits<double>::infinity();
    const Eigenfunction &b = efs.at(t.pair[0]), &c = efs.at(t.pair[1]);
    return stab::threshold_energy(kappa, b.omega, c.omega, b.gamma, c.gamma, t.delta);
}

void write_network (const std::filesystem::path& p, const std::string& header,
                    const std::vector<RadialTriplet>& chosen, const std::vector<RadialTriplet>& net,
                    const std::vector<Key>& by_id, const std::set<Key>& parents,
                    const ModeMap& efs, const stab::KappaCache& kap) {
    std::map<Key, int> id;
    for (size_t i = 0; i < by_id.size(); ++i) id[by_id[i]] = int(i);
    std::ofstream f(p);
    const std::string rule_line = "#" + std::string(98, '-') + "\n";
    char buf[200];
    f << rule_line << header
      << "# selecting legs (parent -> b + c): Delta [c/d], kappa, E_th/E_star\n";
    for (const auto& t : chosen) {
        std::snprintf(buf, sizeof buf, "#   (%d,%+d) -> (%d,%+d) + (%d,%+d)   %+.3e  %+.3e  %.3e\n",
                      t.sum_mode.first, t.sum_mode.second, t.pair[0].first, t.pair[0].second,
                      t.pair[1].first, t.pair[1].second, t.delta / CD, kappa_of(kap, t),
                      e_th(efs, t, kappa_of(kap, t)));
        f << buf;
    }
    f << rule_line << "  2nd_order\n" << rule_line
      << "#  id,   gen,      n,         l,        m,           omega,             gamma,           flin\n"
      << rule_line;
    for (size_t i = 0; i < by_id.size(); ++i) {
        const Eigenfunction& e = efs.at(by_id[i]);
        std::snprintf(buf, sizeof buf, "%5zu %5d %10.1f %9.1f %9.1f %21.10e %17.5e %17.5e\n",
                      i, parents.count(by_id[i]) ? 0 : 1, double(e.n_pg), double(e.l), 0.0,
                      e.omega, e.gamma, 0.0);
        f << buf;
    }
    f << rule_line << "#  id1,   id2,   id3,        kappa\n" << rule_line;
    for (const auto& t : net) {
        const auto k = t.keys();
        std::snprintf(buf, sizeof buf, "%6d %6d %6d %20.10e\n",
                      id.at(k[0]), id.at(k[1]), id.at(k[2]), kappa_of(kap, t));
        f << buf;
    }
    f << rule_line;
}

}  // namespace

int main (int argc, char** argv) {
    CLI::App app{"Parametric networks around one or two parents"};
    std::string model = "models/dsct_M2.0", detail_dir = "detail_wide";
    std::string inlist = "gyre_ad_wide.in", gamma_mode = "tot";
    std::string out = "data/network", cache = "out/four_mode_dsct_M2.0.kappa_cache.tsv";
    std::vector<std::string> nad = {"summary_nad_wide.h5", "summary_nad_wide_hi.h5"};
    std::vector<std::string> parent_s, rules = {"eth"};
    std::vector<double> cuts = {DETUNING_CUT_DIMLESS};
    std::vector<int> ns = {2};
    double closure_dimless = -1.0;
    int l_max = 25, jobs = 4;
    unsigned seed = 1;
    bool closure = true;
    app.add_option("--model", model);
    app.add_option("--detail-dir", detail_dir);
    app.add_option("--inlist", inlist);
    app.add_option("--nad", nad)->expected(-1);
    app.add_option("--gamma", gamma_mode)->check(CLI::IsMember({"rad", "tot"}));
    app.add_option("--parent", parent_s, "l,n of a parent; one or two")->required();
    app.add_option("--rule", rules)->check(CLI::IsMember({"eth", "delta", "dgamma", "random"}));
    app.add_option("-N", ns, "daughter counts");
    app.add_option("--cut", cuts, "candidate |Delta| cuts, units of sqrt(GM/R^3)");
    app.add_option("--closure-cut", closure_dimless, "closure |Delta|; default each --cut");
    app.add_flag("!--no-closure", closure);
    app.add_option("--l-max", l_max);
    app.add_option("--seed", seed, "for --rule random");
    app.add_option("--kappa-cache", cache);
    app.add_option("-j,--jobs", jobs);
    app.add_option("--out", out, "directory");
    CLI11_PARSE(app, argc, argv);
    if (parent_s.size() > 2) { std::fprintf(stderr, "one or two parents\n"); return 1; }

    Model m = load_model(model, detail_dir, inlist, nad);
    ModeMap efs;
    for (auto& [k, e] : m.eigfuncs)
        if (std::isfinite(e.gamma) && k.first <= l_max) efs.emplace(k, e);
    if (gamma_mode == "tot")
        for (auto& [k, e] : efs) e.gamma += gamma_turb(e);
    const double wdyn = m.star->omega_dyn();

    std::set<Key> parents;
    for (const auto& s : parent_s) {
        const Key k = parse_key(s);
        if (!efs.count(k) || efs.at(k).gamma >= 0.0) {
            std::fprintf(stderr, "(%d,%+d) is not a self-excited mode under gamma_%s\n",
                         k.first, k.second, gamma_mode.c_str());
            return 1;
        }
        parents.insert(k);
    }

    // Enumerated once at the widest cut; every narrower one is a subset.
    std::vector<RadialTriplet> all;
    const double cut_max = *std::max_element(cuts.begin(), cuts.end()) * wdyn;
    for (const auto& t : enumerate_triplets(efs, cut_max, l_max, &parents))
        if (efs.at(t.pair[0]).gamma > 0.0 && efs.at(t.pair[1]).gamma > 0.0) all.push_back(t);
    std::printf("%zu modes, %zu parents, %zu candidate pairs at cut %.3g c/d\n",
                efs.size(), parents.size(), all.size(), cut_max / CD);
    if (all.empty()) return 1;

    stab::KappaCache kap = stab::load_kappa_cache(cache);
    auto legs_of = [](const std::vector<RadialTriplet>& ts) {
        std::vector<Leg> out;
        for (const auto& t : ts) { const auto k = t.keys(); out.push_back(stab::leg_of(k[0], k[1], k[2])); }
        return out;
    };
    if (std::count(rules.begin(), rules.end(), "eth"))
        stab::kappa_m000(legs_of(all), efs, cache, kap, jobs);
    std::filesystem::create_directories(out);

    for (double cut_dimless : cuts) {
        std::vector<RadialTriplet> cand;
        for (const auto& t : all) if (std::abs(t.delta) < cut_dimless * wdyn) cand.push_back(t);
        const double closure_cut = (closure_dimless < 0.0 ? cut_dimless : closure_dimless) * wdyn;

        for (const std::string& rule : rules) {
            std::vector<double> score(cand.size());
            std::mt19937 rng(seed);
            std::uniform_real_distribution<double> u;
            for (size_t i = 0; i < cand.size(); ++i) {
                const RadialTriplet& t = cand[i];
                const double gsum = efs.at(t.pair[0]).gamma + efs.at(t.pair[1]).gamma;
                score[i] = rule == "eth"    ? e_th(efs, t, kappa_of(kap, t))
                         : rule == "delta"  ? std::abs(t.delta)
                         : rule == "dgamma" ? std::abs(t.delta) / gsum
                                            : u(rng);
            }
            std::vector<size_t> order(cand.size());
            std::iota(order.begin(), order.end(), size_t(0));
            std::stable_sort(order.begin(), order.end(),
                             [&](size_t a, size_t b) { return score[a] < score[b]; });

            for (int n_daughters : ns) {
                std::vector<Key> daughters;
                std::vector<RadialTriplet> chosen;
                for (size_t i : order) {
                    if (int(daughters.size()) >= n_daughters) break;
                    if (rule == "eth" && !std::isfinite(score[i])) break;
                    chosen.push_back(cand[i]);
                    for (Key k : cand[i].pair)
                        if (std::find(daughters.begin(), daughters.end(), k) == daughters.end())
                            daughters.push_back(k);
                }

                std::vector<RadialTriplet> net = chosen;
                if (closure) {
                    ModeMap sub;
                    for (Key k : parents) sub.emplace(k, efs.at(k));
                    for (Key k : daughters) sub.emplace(k, efs.at(k));
                    net = enumerate_triplets(sub, closure_cut, l_max);
                    const auto have = legs_of(net);
                    const std::set<Leg> in(have.begin(), have.end());
                    for (const auto& t : chosen)         // a closure cut below the cut
                        if (!in.count(legs_of({t})[0])) net.push_back(t);
                }
                stab::kappa_m000(legs_of(net), efs, cache, kap, jobs);
                net.erase(std::remove_if(net.begin(), net.end(), [&](const RadialTriplet& t) {
                              return kappa_of(kap, t) == 0.0; }),
                          net.end());

                std::vector<Key> by_id(parents.begin(), parents.end());
                by_id.insert(by_id.end(), daughters.begin(), daughters.end());
                char name[96], header[400];
                std::snprintf(name, sizeof name, "%s_cut%g_N%d%s.data", rule.c_str(),
                              cut_dimless, n_daughters, closure ? "" : "_open");
                std::snprintf(header, sizeof header,
                              "# network_build %s  gamma_%s  rule %s  N %d  cut %g  closure %s\n",
                              model.c_str(), gamma_mode.c_str(), rule.c_str(), n_daughters,
                              cut_dimless, closure ? std::to_string(closure_cut / wdyn).c_str() : "off");
                const auto p = std::filesystem::path(out) / name;
                write_network(p, header, chosen, net, by_id, parents, efs, kap);
                std::printf("  %-28s %3zu daughters, %3zu legs, %4zu triplets\n",
                            name, daughters.size(), chosen.size(), net.size());
            }
        }
    }
    return 0;
}
