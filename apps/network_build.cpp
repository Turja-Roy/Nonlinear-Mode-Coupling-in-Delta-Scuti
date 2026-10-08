/* Parametric networks around one or two parents, written for network_run:

       ./network_build --parent 6,3 --parent 6,4 --rule eth delta -N 2 4 8 --out data/net
       ./network_run data/net/eth_cut0.15_N8.data --out out/net/eth_cut0.15_N8

   One file per (cut, rule, N) under --out, from one load of the model: the
   load is minutes, a network is milliseconds.

   Candidates are the triplets a parent heads, parent -> b + c with both
   daughters damped and |Delta| < cut; b == c (a -> d + d) included.
   `rule` ranks each parent's pairs, and the parents take turns, best pair
   first, until N daughters are in (one more if the last pair brings two);
   every parent gets at least one pair. An N that brings no new pair (N = 2
   after N = 1) is not written again. Net-growing triplets
   (gamma_b + gamma_c < |gamma_parent|) are kept unless --min-gamma-ratio
   drops them: whether added modes bring their runaway down is the question.

       eth      E_th ascending -- kappa, gamma and Delta together
       eeq      larger daughter energy at the MW25 A7 fixed point ascending;
                infinite for a net-growing triplet. Punishes a weakly damped
                daughter, which must hold ~ |gamma_p|/gamma_d x E_parent there
       delta    |Delta| ascending
       dgamma   |Delta| / (gamma_b + gamma_c) ascending, E_th without kappa
       random   uniform; the control for how much the ranking matters

   --min-gamma-ratio r keeps a pair only if min(gamma_b, gamma_c) >= r |gamma_p|.
   --self-coupled k adds each parent's k best a -> d + d pairs on top of N;
   they are rare (one mode at omega/2, not any pair summing to omega) and
   rank low, so without it they are seldom picked. None exists for odd l_a:
   l_a + 2 l_d must be even.

   Closure then adds every triplet among the chosen modes with |Delta| below
   the closure cut, whichever leg selected them: the equations do not know why
   a mode is in the network. --no-closure keeps the selecting pairs only. */

#include "csv.hpp"
#include "stability.hpp"

#include <CLI/CLI.hpp>

#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>

using stab::TripletKey;

namespace {

Key parse_key (const std::string& s) {
    const size_t c = s.find(',');
    if (c == std::string::npos) throw std::runtime_error("--parent wants l,n: " + s);
    return {std::stoi(s.substr(0, c)), std::stoi(s.substr(c + 1))};
}

double kappa_of (const stab::KappaCache& kap, const RadialTriplet& t) {
    const auto k = t.keys();
    return kap.at(stab::triplet_key(k[0], k[1], k[2])).first;
}

double e_th (const ModeMap& efs, const RadialTriplet& t, double kappa) {
    if (kappa == 0.0) return std::numeric_limits<double>::infinity();
    const Eigenfunction &b = efs.at(t.pair[0]), &c = efs.at(t.pair[1]);
    return stab::threshold_energy(kappa, b.omega, c.omega, b.gamma, c.gamma, t.delta);
}

double e_eq (const ModeMap& efs, const RadialTriplet& t, double kappa) {
    const double inf = std::numeric_limits<double>::infinity();
    const std::array<double, 3> g = {efs.at(t.sum_mode).gamma, efs.at(t.pair[0]).gamma,
                                     efs.at(t.pair[1]).gamma};
    if (kappa == 0.0 || g[0] + g[1] + g[2] <= 0.0) return inf;
    const auto E = stab::equilibrium_energies(kappa, t.omega, g, t.delta);
    const double e = std::max(E[1], E[2]);
    return std::isfinite(e) ? e : inf;
}

double frac_detuning (const RadialTriplet& t) { return t.delta / -t.omega[0]; }

bool self_coupled (const RadialTriplet& t) { return t.pair[0] == t.pair[1]; }

/* Each mode twice, id i at +omega and i + M at -omega, and each triplet as
   (s-, p+, q+) and its mirror (s+, p-, q-): see amplitude.hpp. */
void write_network (const std::filesystem::path& p, const std::string& header,
                    const std::vector<RadialTriplet>& chosen, const std::vector<RadialTriplet>& net,
                    const std::vector<Key>& by_id, const std::set<Key>& parents,
                    const ModeMap& efs, const stab::KappaCache& kap) {
    std::map<Key, int> id;
    const int M = int(by_id.size());
    for (int i = 0; i < M; ++i) id[by_id[size_t(i)]] = i;
    std::ofstream f(p);
    const std::string rule_line = "#" + std::string(98, '-') + "\n";
    char buf[200];
    f << rule_line << header
      << "# selecting pairs (parent -> b + c): Delta [c/d], Delta/omega_parent, kappa, E_th/E_star\n";
    for (const auto& t : chosen) {
        std::snprintf(buf, sizeof buf, "#   (%d,%+d) -> (%d,%+d) + (%d,%+d)   %+.3e  %+.3e  %+.3e  %.3e\n",
                      t.sum_mode.first, t.sum_mode.second, t.pair[0].first, t.pair[0].second,
                      t.pair[1].first, t.pair[1].second, t.delta / CD, frac_detuning(t),
                      kappa_of(kap, t), e_th(efs, t, kappa_of(kap, t)));
        f << buf;
    }
    f << rule_line << "  2nd_order\n" << rule_line
      << "#  id,   gen,      n,         l,        m,           omega,             gamma\n"
      << rule_line;
    for (int s : {+1, -1})
        for (int i = 0; i < M; ++i) {
            const Eigenfunction& e = efs.at(by_id[size_t(i)]);
            std::snprintf(buf, sizeof buf, "%5d %5d %10.1f %9.1f %9.1f %21.10e %17.5e\n",
                          s > 0 ? i : i + M, parents.count(by_id[size_t(i)]) ? 0 : 1,
                          double(e.n_pg), double(e.l), 0.0, s * e.omega, e.gamma);
            f << buf;
        }
    f << rule_line << "#  id1,   id2,   id3,        kappa,       Delta/omega_1\n" << rule_line;
    for (int s : {-1, +1})                       // s: sign of the sum mode's copy
        for (const auto& t : net) {
            const auto k = t.keys();
            const int off_s = s > 0 ? 0 : M, off_p = s > 0 ? M : 0;
            std::snprintf(buf, sizeof buf, "%6d %6d %6d %20.10e %+15.5e\n",
                          id.at(k[0]) + off_s, id.at(k[1]) + off_p, id.at(k[2]) + off_p,
                          kappa_of(kap, t), -s * frac_detuning(t));
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
    double closure_dimless = -1.0, min_gamma_ratio = 0.0;
    int l_max = 25, jobs = 4, n_parents = 1, parent_l_max = 2, n_self = 0;
    unsigned seed = 1;
    bool closure = true;
    app.add_option("--model", model);
    app.add_option("--detail-dir", detail_dir);
    app.add_option("--inlist", inlist);
    app.add_option("--nad", nad)->expected(-1);
    app.add_option("--gamma", gamma_mode)->check(CLI::IsMember({"rad", "tot"}));
    app.add_option("--parent", parent_s, "l,n of a parent; one or two");
    app.add_option("--n-parents", n_parents, "without --parent: the modes of lowest E_th");
    app.add_option("--parent-l-max", parent_l_max, "without --parent: largest l a parent may have");
    app.add_option("--rule", rules)->check(CLI::IsMember({"eth", "eeq", "delta", "dgamma", "random"}));
    app.add_option("--min-gamma-ratio", min_gamma_ratio,
                   "keep pairs with min(gamma_b, gamma_c) >= r |gamma_parent|");
    app.add_option("--self-coupled", n_self, "add each parent's k best a -> d + d pairs");
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
    if (parent_s.size() > 2 || n_parents < 1 || n_parents > 2) {
        std::fprintf(stderr, "one or two parents\n");
        return 1;
    }

    Model m = load_model(model, detail_dir, inlist, nad);
    ModeMap efs;
    for (auto& [k, e] : m.eigfuncs)
        if (std::isfinite(e.gamma) && k.first <= l_max) efs.emplace(k, e);
    if (gamma_mode == "tot")
        for (auto& [k, e] : efs) e.gamma += gamma_turb(e);
    const double wdyn = m.star->omega_dyn();

    stab::KappaCache kap = stab::load_kappa_cache(cache);
    auto keys_of = [](const std::vector<RadialTriplet>& ts) {
        std::vector<TripletKey> out;
        for (const auto& t : ts) { const auto k = t.keys(); out.push_back(stab::triplet_key(k[0], k[1], k[2])); }
        return out;
    };
    // Parent-headed, both daughters damped, at the widest cut; narrower cuts are subsets.
    const double cut_max = *std::max_element(cuts.begin(), cuts.end()) * wdyn;
    auto pairs_of = [&](const std::set<Key>& heads) {
        std::vector<RadialTriplet> out;
        for (const auto& t : enumerate_triplets(efs, cut_max, l_max, &heads)) {
            const double gmin = std::min(efs.at(t.pair[0]).gamma, efs.at(t.pair[1]).gamma);
            if (gmin > 0.0 && gmin >= -min_gamma_ratio * efs.at(t.sum_mode).gamma) out.push_back(t);
        }
        return out;
    };

    /* No --parent: the driven l <= parent_l_max modes easiest to destabilise, i.e. lowest
       min E_th over their damped pairs. Not the most driven: a fast-growing
       high overtone puts its daughters at f/2, inside the driven band, and has
       no damped pair at all. */
    std::vector<RadialTriplet> all;
    if (parent_s.empty()) {
        std::set<Key> driven;
        for (const auto& [k, e] : efs) if (e.gamma < 0.0 && k.first <= parent_l_max) driven.insert(k);
        all = pairs_of(driven);
        stab::kappa_m000(keys_of(all), efs, cache, kap, jobs);
        std::map<Key, std::pair<double, long>> best;           // min E_th, pair count
        for (Key k : driven) best[k] = {std::numeric_limits<double>::infinity(), 0};
        for (const auto& t : all) {
            auto& [e, n] = best[t.sum_mode];
            e = std::min(e, e_th(efs, t, kappa_of(kap, t)));
            ++n;
        }
        std::vector<Key> order(driven.begin(), driven.end());
        std::sort(order.begin(), order.end(),
                  [&](Key a, Key b) { return best[a].first < best[b].first; });
        std::printf("driven l <= %d modes by min E_th over damped pairs:\n", parent_l_max);
        for (Key k : order)
            std::printf("  (%d,%+3d)  f %8.4f c/d  gamma %+.3e  %6ld pairs  min E_th %.3e\n",
                        k.first, k.second, efs.at(k).omega / CD, efs.at(k).gamma,
                        best[k].second, best[k].first);
        for (int i = 0; i < n_parents && i < int(order.size()); ++i)
            if (best[order[size_t(i)]].second)
                parent_s.push_back(std::to_string(order[size_t(i)].first) + ","
                                   + std::to_string(order[size_t(i)].second));
    }
    std::set<Key> parents;
    for (const auto& s : parent_s) {
        const Key k = parse_key(s);
        if (!efs.count(k) || efs.at(k).gamma >= 0.0) {
            std::fprintf(stderr, "(%d,%+d) is not a self-excited mode under gamma_%s\n",
                         k.first, k.second, gamma_mode.c_str());
            return 1;
        }
        parents.insert(k);
        std::printf("parent (%d,%+d)  f %.4f c/d  gamma %+.3e s^-1\n",
                    k.first, k.second, efs.at(k).omega / CD, efs.at(k).gamma);
    }
    if (all.empty()) all = pairs_of(parents);
    all.erase(std::remove_if(all.begin(), all.end(),
                             [&](const RadialTriplet& t) { return !parents.count(t.sum_mode); }),
              all.end());
    std::printf("%zu modes, %zu parents, %zu candidate pairs at cut %.3g c/d\n",
                efs.size(), parents.size(), all.size(), cut_max / CD);
    for (Key k : parents) {
        long n = 0, n_sc = 0;
        double best_sc = std::numeric_limits<double>::infinity();
        for (const auto& t : all) {
            if (t.sum_mode != k) continue;
            ++n;
            if (self_coupled(t)) { ++n_sc; best_sc = std::min(best_sc, std::abs(frac_detuning(t))); }
        }
        std::printf("  (%d,%+d): %ld pairs with min gamma_d >= %g |gamma_p|, %ld of them a -> d + d"
                    " (best |Delta|/omega %.2e)\n", k.first, k.second, n, min_gamma_ratio, n_sc, best_sc);
    }
    if (all.empty()) {
        std::fprintf(stderr, "no damped pair passes: the daughters near f/2 are driven, "
                             "or --min-gamma-ratio leaves none\n");
        return 1;
    }
    if (std::count(rules.begin(), rules.end(), "eth") || std::count(rules.begin(), rules.end(), "eeq"))
        stab::kappa_m000(keys_of(all), efs, cache, kap, jobs);
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
                         : rule == "eeq"    ? e_eq(efs, t, kappa_of(kap, t))
                         : rule == "delta"  ? std::abs(t.delta)
                         : rule == "dgamma" ? std::abs(t.delta) / gsum
                                            : u(rng);
            }
            // Per parent, best first; unusable (infinite score) pairs dropped.
            std::map<Key, std::vector<size_t>> order;
            for (size_t i = 0; i < cand.size(); ++i)
                if (std::isfinite(score[i])) order[cand[i].sum_mode].push_back(i);
            for (auto& [k, o] : order)
                std::stable_sort(o.begin(), o.end(),
                                 [&](size_t a, size_t b) { return score[a] < score[b]; });

            std::vector<int> ns_up(ns);
            std::sort(ns_up.begin(), ns_up.end());
            std::set<TripletKey> prev;                   // pairs of the last file written
            int prev_n = 0;
            for (int n_daughters : ns_up) {
                std::vector<Key> daughters;
                std::vector<RadialTriplet> chosen;
                std::set<size_t> taken;
                auto take = [&](size_t i) {
                    if (!taken.insert(i).second) return;
                    chosen.push_back(cand[i]);
                    for (Key k : cand[i].pair)
                        if (std::find(daughters.begin(), daughters.end(), k) == daughters.end())
                            daughters.push_back(k);
                };
                // Parents take turns; round 0 runs to the end so each has a pair.
                for (size_t r = 0;; ++r) {
                    bool any = false;
                    for (const auto& [k, o] : order) {
                        if (r >= o.size()) continue;
                        if (r > 0 && int(daughters.size()) >= n_daughters) break;
                        take(o[r]);
                        any = true;
                    }
                    if (!any || int(daughters.size()) >= n_daughters) break;
                }
                for (const auto& [k, o] : order) {
                    int n = 0;
                    for (size_t i : o)
                        if (n < n_self && self_coupled(cand[i])) { take(i); ++n; }
                }

                /* Pairs come whole, so a larger N can give the same pairs (N = 1
                   and 2 always do). Write each network once, at its smallest N. */
                const auto ck = keys_of(chosen);
                std::set<TripletKey> now(ck.begin(), ck.end());
                if (now == prev) {
                    std::printf("  N %d: same pairs as N %d, not written\n", n_daughters, prev_n);
                    continue;
                }
                prev = std::move(now);
                prev_n = n_daughters;

                std::vector<RadialTriplet> net = chosen;
                if (closure) {
                    ModeMap sub;
                    for (Key k : parents) sub.emplace(k, efs.at(k));
                    for (Key k : daughters) sub.emplace(k, efs.at(k));
                    net = enumerate_triplets(sub, closure_cut, l_max);
                    const auto have = keys_of(net);
                    const std::set<TripletKey> in(have.begin(), have.end());
                    for (const auto& t : chosen)         // a closure cut below the cut
                        if (!in.count(keys_of({t})[0])) net.push_back(t);
                }
                stab::kappa_m000(keys_of(net), efs, cache, kap, jobs);
                net.erase(std::remove_if(net.begin(), net.end(), [&](const RadialTriplet& t) {
                              return kappa_of(kap, t) == 0.0; }),
                          net.end());

                std::vector<Key> by_id(parents.begin(), parents.end());
                by_id.insert(by_id.end(), daughters.begin(), daughters.end());
                char name[128], header[400], tags[48] = "";
                if (min_gamma_ratio > 0.0)
                    std::snprintf(tags, sizeof tags, "_g%g", min_gamma_ratio);
                if (n_self > 0)
                    std::snprintf(tags + std::strlen(tags), sizeof tags - std::strlen(tags), "_sc%d", n_self);
                std::snprintf(name, sizeof name, "%s_cut%g_N%d%s%s.data", rule.c_str(),
                              cut_dimless, n_daughters, tags, closure ? "" : "_open");
                std::snprintf(header, sizeof header,
                              "# network_build %s  gamma_%s  rule %s  N %d  cut %g  closure %s"
                              "  min_gamma_ratio %g  self_coupled %d\n",
                              model.c_str(), gamma_mode.c_str(), rule.c_str(), n_daughters,
                              cut_dimless, closure ? std::to_string(closure_cut / wdyn).c_str() : "off",
                              min_gamma_ratio, n_self);
                const auto p = std::filesystem::path(out) / name;
                write_network(p, header, chosen, net, by_id, parents, efs, kap);
                std::vector<double> fd;
                for (const auto& t : net) fd.push_back(std::abs(frac_detuning(t)));
                std::sort(fd.begin(), fd.end());
                std::printf("  %-28s %3zu daughters, %3zu pairs, %4zu triplets, "
                            "|Delta|/omega min %.2e median %.2e\n",
                            name, daughters.size(), chosen.size(), net.size(),
                            fd.empty() ? 0.0 : fd.front(), fd.empty() ? 0.0 : fd[fd.size() / 2]);
            }
        }
    }
    return 0;
}
