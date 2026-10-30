/* Networks around one or two parents, written for network_run:

       ./network_build --parent 6,3 -N 2 4 6 8 --out data/net
       ./network_run data/net/cut0.15_d8.data --out out/net/cut0.15_d8

   One file per (cut, N) under --out, from one load of the model: the load is
   minutes, a network is milliseconds.

   Candidates are the triplets a parent heads, parent -> b + c with both
   daughters damped and |Delta| < cut; b == c (a -> d + d) included. Each
   parent's pairs are ranked by E_th, lowest first. Every parent's best pair
   goes in, then the parents take turns, skipping a pair that would bring the
   daughter count past N, until exactly N are in; an N that no set of pairs
   hits is not written. Net-growing triplets (gamma_b + gamma_c < |gamma_parent|)
   are kept unless --min-gamma-ratio drops them: whether added modes bring
   their runaway down is the question.

   --shared (two parents): N counts shared daughters instead, d with a pair
   under each parent, p1 -> d + x and p2 -> d + y, each the best by E_th that
   holds d. Ranked by the larger of the two E_th, lowest first: d couples the
   parents only once both pairs grow. Each d brings both pairs.

   --min-gamma-ratio r keeps a pair only if min(gamma_b, gamma_c) >= r |gamma_p|.
   --self-coupled k adds each parent's k best a -> d + d pairs on top of N;
   they are rare (one mode at omega/2, not any pair summing to omega) and
   rank low, so without it they are seldom picked. None exists for odd l_a:
   l_a + 2 l_d must be even.
   --direct k adds, on top of N, the k best direct legs x + y -> c for each
   pair of parents (x == y: a + a -> c, at 2 omega_a), c damped and
   |Delta| < cut, ranked by the energy c is forced to, E_c = R_c E_x E_y,
   R_c = (s mu)^2 with s = 2 for x != y and 1 for x == y (stab::mu; checked
   against the equations in test_amplitude). A difference leg, c at
   |omega_x - omega_y|, needs the two parents at opposite signs, which signing
   by generation cannot give.

   Closure then adds every triplet among the chosen modes with |Delta| below
   the closure cut, whichever leg selected them: the equations do not know why
   a mode is in the network. network_run signs omega by generation
   (amplitude.hpp), so a triplet resonates only with its sum slot in one
   generation and both pair slots in the other: parametric 0 -> 1 + 1, direct
   0 + 0 -> 1. Any other (a daughter -> daughter + daughter, a difference leg)
   is left out and counted. --no-closure keeps the selecting legs only. */

#include "csv.hpp"
#include "stability.hpp"

#include <CLI/CLI.hpp>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>

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

// Direct leg x + y -> c: E_c / (E_x E_y) with the parents held.
double r_c (const ModeMap& efs, const RadialTriplet& t, double kappa) {
    const Eigenfunction& c = efs.at(t.sum_mode);
    const double s = t.pair[0] == t.pair[1] ? 1.0 : 2.0;
    return std::pow(s * stab::mu(kappa, c.omega, t.delta, c.gamma), 2);
}

double frac_detuning (const RadialTriplet& t) { return t.delta / -t.omega[0]; }

bool self_coupled (const RadialTriplet& t) { return t.pair[0] == t.pair[1]; }

// One row per mode, omega a magnitude: network_run signs it by generation.
void write_network (const std::filesystem::path& p, const std::string& header,
                    const std::vector<RadialTriplet>& chosen, const std::vector<RadialTriplet>& legs,
                    const std::vector<RadialTriplet>& net, const std::vector<Key>& by_id,
                    const std::set<Key>& parents, const ModeMap& efs, const stab::KappaCache& kap) {
    std::map<Key, int> id;
    for (size_t i = 0; i < by_id.size(); ++i) id[by_id[i]] = int(i);
    // Delta as network_run has it, the signed sum: a direct leg's sum slot is the + mode
    auto sign = [&](const RadialTriplet& t) { return parents.count(t.sum_mode) ? 1.0 : -1.0; };
    std::ofstream f(p);
    const std::string rule_line = "#" + std::string(98, '-') + "\n";
    char buf[200];
    f << rule_line << header
      << "# selecting pairs (parent -> b + c): Delta [c/d], Delta/omega_parent, kappa, E_th/E_star"
         " (Weinberg+12 eq. 50), E_eq/E_star (MW25 eq. 6;\n"
      << "#   none: gamma_a + gamma_b + gamma_c <= 0, no stable equilibrium)\n";
    for (const auto& t : chosen) {
        const Eigenfunction &a = efs.at(t.sum_mode), &b = efs.at(t.pair[0]), &c = efs.at(t.pair[1]);
        const double k = kappa_of(kap, t);
        char eq[16] = "none";
        if (a.gamma + b.gamma + c.gamma > 0.0)
            std::snprintf(eq, sizeof eq, "%.3e", stab::equilibrium_energy(k, b.omega, c.omega, a.gamma,
                                                                          b.gamma, c.gamma, t.delta));
        std::snprintf(buf, sizeof buf, "#   (%d,%+d) -> (%d,%+d) + (%d,%+d)   %+.3e  %+.3e  %+.3e  %.3e  %s\n",
                      t.sum_mode.first, t.sum_mode.second, t.pair[0].first, t.pair[0].second,
                      t.pair[1].first, t.pair[1].second, t.delta / CD, frac_detuning(t),
                      k, e_th(efs, t, k), eq);
        f << buf;
    }
    if (!legs.empty())
        f << "# direct legs (x + y -> c): Delta [c/d], Delta/omega_c, kappa, R_c = E_c/(E_x E_y), parents held\n";
    for (const auto& t : legs) {
        const double k = kappa_of(kap, t);
        std::snprintf(buf, sizeof buf, "#   (%d,%+d) + (%d,%+d) -> (%d,%+d)   %+.3e  %+.3e  %+.3e  %.3e\n",
                      t.pair[0].first, t.pair[0].second, t.pair[1].first, t.pair[1].second,
                      t.sum_mode.first, t.sum_mode.second, sign(t) * t.delta / CD,
                      sign(t) * frac_detuning(t), k, r_c(efs, t, k));
        f << buf;
    }
    f << rule_line << "  2nd_order\n" << rule_line
      << "#  id,   gen,      n,         l,        m,           omega,             gamma\n"
      << rule_line;
    for (size_t i = 0; i < by_id.size(); ++i) {
        const Eigenfunction& e = efs.at(by_id[i]);
        std::snprintf(buf, sizeof buf, "%5zu %5d %10.1f %9.1f %9.1f %21.10e %17.5e\n",
                      i, parents.count(by_id[i]) ? 0 : 1, double(e.n_pg), double(e.l), 0.0,
                      e.omega, e.gamma);
        f << buf;
    }
    f << rule_line << "#  id1,   id2,   id3,        kappa,       Delta/|omega_1|\n" << rule_line;
    for (const auto& t : net) {
        const auto k = t.keys();
        std::snprintf(buf, sizeof buf, "%6d %6d %6d %20.10e %+15.5e\n",
                      id.at(k[0]), id.at(k[1]), id.at(k[2]), kappa_of(kap, t), sign(t) * frac_detuning(t));
        f << buf;
    }
    f << rule_line;
}

}  // namespace

int main (int argc, char** argv) {
    CLI::App app{"Networks around one or two parents"};
    std::string model = "models/dsct_M2.0", detail_dir = "detail_wide";
    std::string inlist = "gyre_ad_wide.in", gamma_mode = "tot";
    std::string out = "data/network", cache = "out/four_mode_dsct_M2.0.kappa_cache.tsv";
    std::vector<std::string> nad = {"summary_nad_wide.h5", "summary_nad_wide_hi.h5"};
    std::vector<std::string> parent_s;
    std::vector<double> cuts = {DETUNING_CUT_DIMLESS};
    std::vector<int> ns = {2};
    double closure_dimless = -1.0, min_gamma_ratio = 0.0;
    int l_max = 25, jobs = 4, n_parents = 1, parent_l_max = 2, n_self = 0, n_direct = 0;
    bool closure = true, shared = false;
    app.add_option("--model", model);
    app.add_option("--detail-dir", detail_dir);
    app.add_option("--inlist", inlist);
    app.add_option("--nad", nad)->expected(-1);
    app.add_option("--gamma", gamma_mode)->check(CLI::IsMember({"rad", "tot"}));
    app.add_option("--parent", parent_s, "l,n of a parent; one or two");
    app.add_option("--n-parents", n_parents, "without --parent: the modes of lowest E_th");
    app.add_option("--parent-l-max", parent_l_max, "without --parent: largest l a parent may have");
    app.add_option("--min-gamma-ratio", min_gamma_ratio,
                   "keep pairs with min(gamma_b, gamma_c) >= r |gamma_parent|");
    app.add_option("--self-coupled", n_self, "add each parent's k best a -> d + d pairs");
    app.add_option("--direct", n_direct, "add the k best direct legs x + y -> c per pair of parents");
    app.add_flag("--shared", shared, "two parents: N counts daughters shared by both, each with its two pairs");
    app.add_option("-N", ns, "daughter counts, exact (shared daughters with --shared)");
    app.add_option("--cut", cuts, "candidate |Delta| cuts, units of sqrt(GM/R^3)");
    app.add_option("--closure-cut", closure_dimless, "closure |Delta|; default each --cut");
    app.add_flag("!--no-closure", closure);
    app.add_option("--l-max", l_max);
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
    if (shared && parents.size() != 2) {
        std::fprintf(stderr, "--shared needs two parents\n");
        return 1;
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
    stab::kappa_m000(keys_of(all), efs, cache, kap, jobs);

    // Direct legs at the widest cut: damped c near omega_x + omega_y, both pair slots parents.
    std::vector<RadialTriplet> legs_all;
    if (n_direct > 0) {
        std::set<Key> near;
        for (Key x : parents)
            for (Key y : parents)
                for (const auto& [k, e] : efs)
                    if (e.gamma > 0.0 && std::abs(e.omega - efs.at(x).omega - efs.at(y).omega) < cut_max)
                        near.insert(k);
        for (const auto& t : enumerate_triplets(efs, cut_max, l_max, &near))
            if (parents.count(t.pair[0]) && parents.count(t.pair[1])) legs_all.push_back(t);
        stab::kappa_m000(keys_of(legs_all), efs, cache, kap, jobs);
        std::printf("%zu direct legs at cut %.3g c/d\n", legs_all.size(), cut_max / CD);
    }
    std::filesystem::create_directories(out);

    for (double cut_dimless : cuts) {
        std::vector<RadialTriplet> cand;
        for (const auto& t : all) if (std::abs(t.delta) < cut_dimless * wdyn) cand.push_back(t);
        const double closure_cut = (closure_dimless < 0.0 ? cut_dimless : closure_dimless) * wdyn;

        std::vector<double> score(cand.size());
        for (size_t i = 0; i < cand.size(); ++i) score[i] = e_th(efs, cand[i], kappa_of(kap, cand[i]));
        // Per parent, best first; unusable (infinite score) pairs dropped.
        std::map<Key, std::vector<size_t>> order;
        for (size_t i = 0; i < cand.size(); ++i)
            if (std::isfinite(score[i])) order[cand[i].sum_mode].push_back(i);
        for (auto& [k, o] : order)
            std::stable_sort(o.begin(), o.end(), [&](size_t a, size_t b) { return score[a] < score[b]; });

        // Each parent's best pair, then the parents in turn.
        std::vector<size_t> turns;
        for (size_t r = 0;; ++r) {
            const size_t before = turns.size();
            for (const auto& [k, o] : order) if (r < o.size()) turns.push_back(o[r]);
            if (turns.size() == before) break;
        }

        // Shared daughters: under each parent the best pair holding d, scored by the larger E_th.
        std::vector<std::pair<double, std::array<size_t, 2>>> ranked;
        if (shared) {
            std::map<Key, std::map<Key, size_t>> best;
            for (const auto& [p, o] : order)
                for (size_t i : o)
                    for (Key d : cand[i].pair) best[d].emplace(p, i);
            for (const auto& [d, by_p] : best)
                if (by_p.size() == 2) {
                    const size_t i = by_p.begin()->second, j = by_p.rbegin()->second;
                    ranked.push_back({std::max(score[i], score[j]), {i, j}});
                }
            std::sort(ranked.begin(), ranked.end());
            std::printf("  cut %g: %zu shared daughters, best larger-E_th %.3e\n", cut_dimless, ranked.size(),
                        ranked.empty() ? std::nan("") : ranked[0].first);
        }

        // The k best direct legs per pair of parents, by R_c.
        std::map<std::pair<Key, Key>, std::vector<RadialTriplet>> by_parents;
        for (const auto& t : legs_all)
            if (std::abs(t.delta) < cut_dimless * wdyn) by_parents[{t.pair[0], t.pair[1]}].push_back(t);
        std::vector<RadialTriplet> legs;
        for (auto& [xy, ts] : by_parents) {
            std::sort(ts.begin(), ts.end(), [&](const RadialTriplet& a, const RadialTriplet& b) {
                return r_c(efs, a, kappa_of(kap, a)) > r_c(efs, b, kappa_of(kap, b)); });
            legs.insert(legs.end(), ts.begin(), ts.begin() + std::min(ts.size(), size_t(n_direct)));
        }

        std::vector<int> ns_up(ns);
        std::sort(ns_up.begin(), ns_up.end());
        for (int n_daughters : ns_up) {
            std::vector<Key> daughters;
            std::vector<RadialTriplet> chosen;
            std::set<size_t> taken;
            auto fresh = [&](size_t i) {
                const std::set<Key> d(cand[i].pair.begin(), cand[i].pair.end());
                return int(std::count_if(d.begin(), d.end(), [&](Key k) {
                    return std::find(daughters.begin(), daughters.end(), k) == daughters.end(); }));
            };
            auto take = [&](size_t i) {
                if (!taken.insert(i).second) return;
                chosen.push_back(cand[i]);
                for (Key k : cand[i].pair)
                    if (std::find(daughters.begin(), daughters.end(), k) == daughters.end())
                        daughters.push_back(k);
            };
            if (shared) {
                if (n_daughters > int(ranked.size())) {
                    std::printf("  %d shared: only %zu within the cut, not written\n", n_daughters, ranked.size());
                    continue;
                }
                for (int j = 0; j < n_daughters; ++j)
                    for (size_t i : ranked[size_t(j)].second) take(i);
            } else {
                for (size_t j = 0; j < turns.size(); ++j)
                    if (j < order.size() || (int(daughters.size()) < n_daughters
                                             && int(daughters.size()) + fresh(turns[j]) <= n_daughters))
                        take(turns[j]);
                if (int(daughters.size()) != n_daughters) {
                    std::printf("  N %d: no exact fit within the cut (%zu daughters), not written\n",
                                n_daughters, daughters.size());
                    continue;
                }
            }
            for (const auto& [k, o] : order) {
                int n = 0;
                for (size_t i : o)
                    if (n < n_self && self_coupled(cand[i])) { take(i); ++n; }
            }
            for (const auto& t : legs)
                if (std::find(daughters.begin(), daughters.end(), t.sum_mode) == daughters.end())
                    daughters.push_back(t.sum_mode);

            std::vector<RadialTriplet> selecting = chosen;
            selecting.insert(selecting.end(), legs.begin(), legs.end());
            std::vector<RadialTriplet> net = selecting;
            if (closure) {
                ModeMap sub;
                for (Key k : parents) sub.emplace(k, efs.at(k));
                for (Key k : daughters) sub.emplace(k, efs.at(k));
                net = enumerate_triplets(sub, closure_cut, l_max);
                const auto have = keys_of(net);
                const std::set<TripletKey> in(have.begin(), have.end());
                for (const auto& t : selecting)          // a closure cut below the cut
                    if (!in.count(keys_of({t})[0])) net.push_back(t);
            }
            auto gen = [&](Key k) { return parents.count(k) ? 0 : 1; };
            const auto off = std::remove_if(net.begin(), net.end(), [&](const RadialTriplet& t) {
                return gen(t.pair[0]) != gen(t.pair[1]) || gen(t.sum_mode) == gen(t.pair[0]); });
            const long n_off = long(net.end() - off);
            net.erase(off, net.end());
            stab::kappa_m000(keys_of(net), efs, cache, kap, jobs);
            net.erase(std::remove_if(net.begin(), net.end(), [&](const RadialTriplet& t) {
                          return kappa_of(kap, t) == 0.0; }),
                      net.end());

            std::vector<Key> by_id(parents.begin(), parents.end());
            by_id.insert(by_id.end(), daughters.begin(), daughters.end());
            char name[128], header[400], tags[48] = "", lead[24] = "";
            if (min_gamma_ratio > 0.0)
                std::snprintf(tags, sizeof tags, "_g%g", min_gamma_ratio);
            if (n_self > 0)
                std::snprintf(tags + std::strlen(tags), sizeof tags - std::strlen(tags), "_sc%d", n_self);
            if (n_direct > 0)
                std::snprintf(tags + std::strlen(tags), sizeof tags - std::strlen(tags), "_dir%d", n_direct);
            if (shared) std::snprintf(lead, sizeof lead, "shared%d_", n_daughters);
            std::snprintf(name, sizeof name, "cut%g_%sd%zu%s%s.data", cut_dimless, lead, daughters.size(),
                          tags, closure ? "" : "_open");
            std::snprintf(header, sizeof header,
                          "# network_build %s  gamma_%s  daughters %zu  cut %g  closure %s"
                          "  min_gamma_ratio %g  self_coupled %d  direct %d  shared %d\n",
                          model.c_str(), gamma_mode.c_str(), daughters.size(), cut_dimless,
                          closure ? std::to_string(closure_cut / wdyn).c_str() : "off",
                          min_gamma_ratio, n_self, n_direct, shared ? n_daughters : 0);
            const auto p = std::filesystem::path(out) / name;
            write_network(p, header, chosen, legs, net, by_id, parents, efs, kap);
            std::vector<double> fd;
            for (const auto& t : net) fd.push_back(std::abs(frac_detuning(t)));
            std::sort(fd.begin(), fd.end());
            std::printf("  %-36s %3zu daughters, %3zu pairs, %zu direct, %4zu triplets (%ld off-generation "
                        "dropped), |Delta|/omega min %.2e median %.2e\n",
                        name, daughters.size(), chosen.size(), legs.size(), net.size(), n_off,
                        fd.empty() ? 0.0 : fd.front(), fd.empty() ? 0.0 : fd[fd.size() / 2]);
        }
    }
    return 0;
}
