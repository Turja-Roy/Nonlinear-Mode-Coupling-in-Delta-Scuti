/* Network ladders around one or two parents, written for network_run:

       ./network_build --parent 2,-2 --out data/net
       ./network_run data/net/l2_n-2/d8.data --out out/net/l2_n-2/d8

   One model load for every ladder: the load is minutes, a network is milliseconds.

   Candidate triplets of a set of parents, all with |Delta| < cut and kappa != 0:
     pair       p -> b + c, b and c damped (b == c included), ranked by E_th (Weinberg+12 eq. 50);
     harmonic   x + x -> c, c damped, ranked by E_eq of the triplet (the hold energy E_h);
     sum        x + y -> c (two parents), ranked by E_eq = E_x + E_y at its fixed point.
   Lowest first, each parent's lowest one ahead of the rest. A ladder takes them in this order; a
   triplet brings 0, 1 or 2 new daughters, and d<N>.data holds every triplet taken before the one that
   brings daughter N + 1, for each N up to --n-max. Each network contains the one before.

   --parent l,n (one or two): one ladder, into <out>/<parents>/.
   --each-parent: one ladder per driven l <= --parent-l-max mode, into <out>/l<l>_n<n>/.
   --pairs k: the k pairs of parents, each held alone by one of its triplets (a pair with an E_eq,
     or a harmonic triplet with |gamma_x| / gamma_c < 1/2), with the lowest E_th or E_eq of a
     triplet coupling them (a shared daughter: the larger E_th of its two pairs; a sum triplet: E_eq).
   --only-eeq keeps only pairs with gamma_p + gamma_b + gamma_c > 0 (an E_eq exists); a ladder that
   it leaves unchanged is not written.
   --min-gamma-ratio r keeps a pair only if min(gamma_b, gamma_c) >= r |gamma_p|.
   --self-coupled k puts each parent's k lowest a -> d + d pairs ahead of the rest.
   A difference triplet, c at |omega_x - omega_y|, needs the two parents at opposite signs, which
   signing by generation cannot give.

   --list-direct: every x + y -> c with x, y driven (l <= --parent-l-max), c damped and |Delta| < cut,
   ranked by E_eq, then by mu, with each parent's damped pairs and how many have an E_eq. No network.

   Closure then adds every triplet among the chosen modes with |Delta| below the closure cut: the
   equations do not know why a mode is in the network. network_run signs omega by generation
   (amplitude.hpp), so only parametric 0 -> 1 + 1 and direct 0 + 0 -> 1 triplets resonate; others
   are left out and counted. --no-closure keeps the selecting triplets only. */

#include "csv.hpp"
#include "stability.hpp"

#include <CLI/CLI.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>

using stab::TripletKey;

namespace {

constexpr double INF = std::numeric_limits<double>::infinity();

Key parse_key (const std::string& s) {
    const size_t c = s.find(',');
    if (c == std::string::npos) throw std::runtime_error("--parent wants l,n: " + s);
    return {std::stoi(s.substr(0, c)), std::stoi(s.substr(c + 1))};
}

std::string name_of (Key k) {
    char b[32];
    std::snprintf(b, sizeof b, "l%d_n%+d", k.first, k.second);
    return b;
}

double kappa_of (const stab::KappaCache& kap, const RadialTriplet& t) {
    const auto k = t.keys();
    return kap.at(stab::triplet_key(k[0], k[1], k[2])).first;
}

double e_th (const ModeMap& efs, const RadialTriplet& t, double kappa) {
    if (kappa == 0.0) return INF;
    const Eigenfunction &b = efs.at(t.pair[0]), &c = efs.at(t.pair[1]);
    return stab::threshold_energy(kappa, b.omega, c.omega, b.gamma, c.gamma, t.delta);
}

// x + y -> c: MW23 eq. 9
double mu_direct (const ModeMap& efs, const RadialTriplet& t, double kappa) {
    const Eigenfunction& c = efs.at(t.sum_mode);
    return stab::mu(kappa, c.omega, t.delta, c.gamma);
}

// x + y -> c: E_x + E_y at the fixed point; x == y: E_x, the hold energy
double e_eq_direct (const ModeMap& efs, const RadialTriplet& t, double kappa) {
    const Eigenfunction &c = efs.at(t.sum_mode), &x = efs.at(t.pair[0]), &y = efs.at(t.pair[1]);
    const bool self = t.pair[0] == t.pair[1];
    const auto E = stab::equilibrium_energies(kappa, {c.omega, -x.omega, -y.omega}, {c.gamma, x.gamma, y.gamma},
                                              t.delta, 1.0, self);
    return self ? E[1] : E[1] + E[2];
}

double frac_detuning (const RadialTriplet& t) { return t.delta / -t.omega[0]; }

bool self_coupled (const RadialTriplet& t) { return t.pair[0] == t.pair[1]; }

bool has_eeq (const ModeMap& efs, const RadialTriplet& t) {
    return efs.at(t.sum_mode).gamma + efs.at(t.pair[0]).gamma + efs.at(t.pair[1]).gamma > 0.0;
}

// held alone: a pair with an E_eq, or a harmonic triplet with |gamma_x| / gamma_c < 1/2
bool holds (const ModeMap& efs, const RadialTriplet& t, bool direct) {
    if (!direct) return has_eeq(efs, t);
    return self_coupled(t) && -efs.at(t.pair[0]).gamma < 0.5 * efs.at(t.sum_mode).gamma;
}

bool involves (const RadialTriplet& t, Key p) { return t.sum_mode == p || t.pair[0] == p || t.pair[1] == p; }

struct Step { std::vector<Key> daughters; std::vector<RadialTriplet> chosen; };

std::vector<Step> ladder (const std::vector<RadialTriplet>& ranked, const std::set<Key>& parents, int n_max) {
    std::vector<Step> steps;
    Step s;
    for (const auto& t : ranked) {
        std::vector<Key> fresh;
        for (Key k : t.keys())
            if (!parents.count(k) && std::find(s.daughters.begin(), s.daughters.end(), k) == s.daughters.end()
                && std::find(fresh.begin(), fresh.end(), k) == fresh.end())
                fresh.push_back(k);
        if (!fresh.empty() && !s.chosen.empty()) steps.push_back(s);
        if (int(s.daughters.size() + fresh.size()) > n_max) return steps;
        s.daughters.insert(s.daughters.end(), fresh.begin(), fresh.end());
        s.chosen.push_back(t);
    }
    if (!s.chosen.empty()) steps.push_back(s);
    return steps;
}

std::string triplet_line (const ModeMap& efs, const stab::KappaCache& kap, const std::set<Key>& parents,
                          const RadialTriplet& t) {
    const double k = kappa_of(kap, t);
    char buf[200];
    if (parents.count(t.sum_mode)) {
        char eq[16] = "none";
        const Eigenfunction &a = efs.at(t.sum_mode), &b = efs.at(t.pair[0]), &c = efs.at(t.pair[1]);
        if (has_eeq(efs, t))
            std::snprintf(eq, sizeof eq, "%.3e", stab::equilibrium_energy(k, b.omega, c.omega, a.gamma,
                                                                          b.gamma, c.gamma, t.delta));
        std::snprintf(buf, sizeof buf, "#   pair      (%d,%+d) -> (%d,%+d) + (%d,%+d)   %+.3e  %+.3e  %+.3e  %.3e  %s\n",
                      t.sum_mode.first, t.sum_mode.second, t.pair[0].first, t.pair[0].second,
                      t.pair[1].first, t.pair[1].second, t.delta / CD, frac_detuning(t), k, e_th(efs, t, k), eq);
    } else {
        std::snprintf(buf, sizeof buf, "#   %-9s (%d,%+d) + (%d,%+d) -> (%d,%+d)   %+.3e  %+.3e  %+.3e  %.3e  mu %.3e\n",
                      self_coupled(t) ? "harmonic" : "sum", t.pair[0].first, t.pair[0].second, t.pair[1].first,
                      t.pair[1].second, t.sum_mode.first, t.sum_mode.second, -t.delta / CD, -frac_detuning(t), k,
                      e_eq_direct(efs, t, k), mu_direct(efs, t, k));
    }
    return buf;
}

// One row per mode, omega a magnitude: network_run signs it by generation.
void write_network (const std::filesystem::path& p, const std::string& header, const Step& s,
                    const std::vector<RadialTriplet>& net, const std::vector<Key>& by_id,
                    const std::set<Key>& parents, const ModeMap& efs, const stab::KappaCache& kap) {
    std::map<Key, int> id;
    for (size_t i=0 ; i<by_id.size() ; i++) id[by_id[i]] = int(i);
    // Delta as network_run has it, the signed sum: a direct triplet's sum slot is the + mode
    auto sign = [&](const RadialTriplet& t) { return parents.count(t.sum_mode) ? 1.0 : -1.0; };
    std::map<Key, std::set<Key>> under;
    for (const auto& t : s.chosen)
        if (parents.count(t.sum_mode)) for (Key d : t.pair) under[d].insert(t.sum_mode);
    std::ofstream f(p);
    const std::string rule_line = "#" + std::string(98, '-') + "\n";
    char buf[200];
    f << rule_line << header
      << "# selecting triplets, in the order they entered: Delta [c/d], Delta/omega, kappa, then\n"
      << "#   pair: E_th/E_star (Weinberg+12 eq. 50), E_eq/E_star (MW25 eq. 6; none: gamma_a + gamma_b + gamma_c <= 0)\n"
      << "#   harmonic, sum: E_eq/E_star at the triplet's fixed point, mu (MW23 eq. 9)\n";
    for (const auto& t : s.chosen) f << triplet_line(efs, kap, parents, t);
    f << "# shared daughters:";
    for (const auto& [d, ps] : under)
        if (ps.size() > 1) { std::snprintf(buf, sizeof buf, " (%d,%+d)", d.first, d.second); f << buf; }
    f << "\n" << rule_line << "  2nd_order\n" << rule_line
      << "#  id,   gen,      n,         l,        m,           omega,             gamma\n"
      << rule_line;
    for (size_t i=0 ; i<by_id.size() ; i++) {
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
    CLI::App app{"Network ladders around one or two parents"};
    std::string model = "models/dsct_M2.0", detail_dir = "detail_wide";
    std::string inlist = "gyre_ad_wide.in", gamma_mode = "tot";
    std::string out = "data/network", cache = "out/four_mode_dsct_M2.0.kappa_cache.tsv";
    std::vector<std::string> nad = {"summary_nad_wide.h5", "summary_nad_wide_hi.h5"};
    std::vector<std::string> parent_s;
    double cut_dimless = DETUNING_CUT_DIMLESS, closure_dimless = -1.0, min_gamma_ratio = 0.0;
    int l_max = 25, jobs = 4, parent_l_max = 2, n_self = 0, n_pairs = 0, n_max = 24;
    bool closure = true, each_parent = false, only_eeq = false, list_direct = false;
    app.add_option("--model", model);
    app.add_option("--detail-dir", detail_dir);
    app.add_option("--inlist", inlist);
    app.add_option("--nad", nad)->expected(-1);
    app.add_option("--gamma", gamma_mode)->check(CLI::IsMember({"rad", "tot"}));
    app.add_option("--parent", parent_s, "l,n of a parent; one or two");
    app.add_flag("--each-parent", each_parent, "one ladder per driven l <= --parent-l-max mode");
    app.add_option("--pairs", n_pairs, "ladders for k pairs of parents, each held alone, lowest E_th or E_eq of a triplet coupling them first");
    app.add_option("--parent-l-max", parent_l_max, "largest l of a parent for --each-parent, --pairs, --list-direct");
    app.add_option("--n-max", n_max, "largest daughter count, default 24");
    app.add_flag("--only-eeq", only_eeq, "only pairs with gamma_p + gamma_b + gamma_c > 0");
    app.add_option("--min-gamma-ratio", min_gamma_ratio,
                   "keep pairs with min(gamma_b, gamma_c) >= r |gamma_parent|");
    app.add_option("--self-coupled", n_self, "put each parent's k lowest a -> d + d pairs first");
    app.add_flag("--list-direct", list_direct, "rank x + y -> c over all driven x, y; no network");
    app.add_option("--cut", cut_dimless, "candidate |Delta| cut, units of sqrt(GM/R^3)");
    app.add_option("--closure-cut", closure_dimless, "closure |Delta|; default --cut");
    app.add_flag("!--no-closure", closure);
    app.add_option("--l-max", l_max);
    app.add_option("--kappa-cache", cache);
    app.add_option("-j,--jobs", jobs);
    app.add_option("--out", out, "directory");
    CLI11_PARSE(app, argc, argv);
    if (parent_s.size() > 2 || (parent_s.empty() && !each_parent && !n_pairs && !list_direct)) {
        std::fprintf(stderr, "one or two --parent, or --each-parent, --pairs, --list-direct\n");
        return 1;
    }

    Model m = load_model(model, detail_dir, inlist, nad);
    ModeMap efs;
    for (auto& [k, e] : m.eigfuncs)
        if (std::isfinite(e.gamma) && k.first <= l_max) efs.emplace(k, e);
    if (gamma_mode == "tot")
        for (auto& [k, e] : efs) e.gamma += gamma_turb(e);
    const double wdyn = m.star->omega_dyn(), cut = cut_dimless * wdyn;
    const double closure_cut = (closure_dimless < 0.0 ? cut_dimless : closure_dimless) * wdyn;

    stab::KappaCache kap = stab::load_kappa_cache(cache);
    auto keys_of = [](const std::vector<RadialTriplet>& ts) {
        std::vector<TripletKey> out;
        for (const auto& t : ts) { const auto k = t.keys(); out.push_back(stab::triplet_key(k[0], k[1], k[2])); }
        return out;
    };

    std::set<Key> driven;
    for (const auto& [k, e] : efs) if (e.gamma < 0.0 && k.first <= parent_l_max) driven.insert(k);
    for (const auto& s : parent_s) {
        const Key k = parse_key(s);
        if (!efs.count(k) || efs.at(k).gamma >= 0.0) {
            std::fprintf(stderr, "(%d,%+d) is not a self-excited mode under gamma_%s\n",
                         k.first, k.second, gamma_mode.c_str());
            return 1;
        }
        driven.insert(k);
    }

    // Parent-headed pairs, both daughters damped.
    std::vector<RadialTriplet> pairs;
    for (const auto& t : enumerate_triplets(efs, cut, l_max, &driven)) {
        const double gmin = std::min(efs.at(t.pair[0]).gamma, efs.at(t.pair[1]).gamma);
        if (gmin > 0.0 && gmin >= -min_gamma_ratio * efs.at(t.sum_mode).gamma) pairs.push_back(t);
    }
    // Direct: damped c near omega_x + omega_y, both pair slots driven.
    std::set<Key> near;
    for (Key x : driven)
        for (Key y : driven)
            for (const auto& [k, e] : efs)
                if (e.gamma > 0.0 && std::abs(e.omega - efs.at(x).omega - efs.at(y).omega) < cut) near.insert(k);
    std::vector<RadialTriplet> direct;
    for (const auto& t : enumerate_triplets(efs, cut, l_max, &near))
        if (driven.count(t.pair[0]) && driven.count(t.pair[1])) direct.push_back(t);
    stab::kappa_m000(keys_of(pairs), efs, cache, kap, jobs);
    stab::kappa_m000(keys_of(direct), efs, cache, kap, jobs);
    auto nonzero = [&](std::vector<RadialTriplet>& ts) {
        ts.erase(std::remove_if(ts.begin(), ts.end(), [&](const RadialTriplet& t) { return kappa_of(kap, t) == 0.0; }),
                 ts.end());
    };
    nonzero(pairs);
    nonzero(direct);
    std::printf("%zu modes, %zu driven, %zu pairs, %zu direct triplets at cut %.3g c/d\n",
                efs.size(), driven.size(), pairs.size(), direct.size(), cut / CD);

    if (list_direct) {
        std::map<Key, std::array<long, 2>> np;          // damped pairs, those with E_eq
        for (const auto& t : pairs) {
            auto& n = np[t.sum_mode];
            ++n[0];
            if (has_eeq(efs, t)) ++n[1];
        }
        std::vector<size_t> order(direct.size());
        std::vector<double> eeq(direct.size()), mu(direct.size());
        for (size_t i=0 ; i<direct.size() ; i++) {
            order[i] = i;
            eeq[i] = e_eq_direct(efs, direct[i], kappa_of(kap, direct[i]));
            mu[i] = mu_direct(efs, direct[i], kappa_of(kap, direct[i]));
        }
        std::printf("E_eq: E_x + E_y at the fixed point (x == y: E_x) [E_star]; mu: MW23 eq. 9\n");
        for (const bool by_eeq : {true, false}) {
            std::stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) {
                return by_eeq ? eeq[i] < eeq[j] : mu[i] > mu[j]; });
            std::printf("\n%s\n#  x        y        c          f_x      f_y      f_c   Delta[c/d]  gamma_x     gamma_y"
                        "     gamma_c     kappa       E_eq       mu         pairs_x eq_x  pairs_y eq_y\n",
                        by_eeq ? "ranked by E_eq, lowest first" : "ranked by mu, highest first");
            for (size_t i : order) {
                const auto& t = direct[i];
                const Eigenfunction &x = efs.at(t.pair[0]), &y = efs.at(t.pair[1]), &c = efs.at(t.sum_mode);
                const auto &nx = np[t.pair[0]], &ny = np[t.pair[1]];
                std::printf("  (%d,%+d)  (%d,%+d)  (%d,%+d)  %8.3f %8.3f %8.3f  %+.3e  %+.3e  %+.3e  %+.3e  %+.3e"
                            "  %.3e  %.3e  %6ld %5ld  %6ld %5ld\n",
                            t.pair[0].first, t.pair[0].second, t.pair[1].first, t.pair[1].second,
                            t.sum_mode.first, t.sum_mode.second, x.omega / CD, y.omega / CD, c.omega / CD,
                            t.delta / CD, x.gamma, y.gamma, c.gamma, kappa_of(kap, t), eeq[i], mu[i],
                            nx[0], nx[1], ny[0], ny[1]);
            }
        }
        return 0;
    }

    // Ranking energy of a candidate: E_th of a pair, E_eq of a direct triplet.
    std::map<TripletKey, double> rank_e;
    for (const auto& t : pairs) rank_e[keys_of({t})[0]] = e_th(efs, t, kappa_of(kap, t));
    for (const auto& t : direct) rank_e[keys_of({t})[0]] = e_eq_direct(efs, t, kappa_of(kap, t));
    auto energy = [&](const RadialTriplet& t) { return rank_e.at(keys_of({t})[0]); };

    auto candidates = [&](const std::set<Key>& ps, bool eq_only) {
        std::vector<RadialTriplet> ts;
        for (const auto& t : pairs)
            if (ps.count(t.sum_mode) && (!eq_only || has_eeq(efs, t))) ts.push_back(t);
        for (const auto& t : direct)
            if (ps.count(t.pair[0]) && ps.count(t.pair[1])) ts.push_back(t);
        std::stable_sort(ts.begin(), ts.end(), [&](const RadialTriplet& a, const RadialTriplet& b) {
            return energy(a) < energy(b); });
        std::vector<size_t> first;
        for (Key p : ps) {
            for (size_t i=0 ; i<ts.size() ; i++) if (involves(ts[i], p)) { first.push_back(i); break; }
            int n = 0;
            for (size_t i=0 ; i<ts.size() && n<n_self ; i++)
                if (ts[i].sum_mode == p && self_coupled(ts[i])) { first.push_back(i); n++; }
        }
        std::sort(first.begin(), first.end());
        first.erase(std::unique(first.begin(), first.end()), first.end());
        std::vector<RadialTriplet> ranked;
        for (size_t i : first) ranked.push_back(ts[i]);
        for (size_t i=0 ; i<ts.size() ; i++)
            if (!std::binary_search(first.begin(), first.end(), i)) ranked.push_back(ts[i]);
        return ranked;
    };

    // Every driven mode: what can hold it.
    std::map<Key, bool> held;
    std::map<Key, std::map<Key, double>> best;   // parent -> daughter -> lowest E_th of its pairs holding it
    std::printf("\ndriven l <= %d modes: damped pairs, those with E_eq, harmonic triplets with |gamma_x|/gamma_c < 1/2\n"
                "#  mode      f [c/d]   gamma [1/s]   pairs   E_eq  harmonic   lowest E_th  lowest E_h   held by\n",
                parent_l_max);
    for (Key k : driven) {
        long n = 0, n_eq = 0, n_h = 0;
        double eth = INF, eh = INF;
        for (const auto& t : pairs)
            if (t.sum_mode == k) {
                ++n;
                if (has_eeq(efs, t)) ++n_eq;
                eth = std::min(eth, energy(t));
                for (Key d : t.pair) {
                    auto it = best[k].emplace(d, INF).first;
                    it->second = std::min(it->second, energy(t));
                }
            }
        for (const auto& t : direct)
            if (t.pair[0] == k && t.pair[1] == k && holds(efs, t, true)) { ++n_h; eh = std::min(eh, energy(t)); }
        held[k] = n_eq || n_h;
        std::printf("  (%d,%+d)  %8.4f   %+.3e   %5ld  %5ld  %8ld   %.3e   %.3e   %s\n", k.first, k.second,
                    efs.at(k).omega / CD, efs.at(k).gamma, n, n_eq, n_h, eth, eh,
                    n_eq ? "a pair" : n_h ? "a harmonic triplet" : n ? "none (pairs without E_eq)"
                                                                     : "none (no triplet within the cut)");
    }

    std::vector<std::vector<Key>> sets;
    if (!parent_s.empty()) {
        sets.emplace_back();
        for (const auto& s : parent_s) sets.back().push_back(parse_key(s));
    }
    if (each_parent) for (Key k : driven) if (k.first <= parent_l_max) sets.push_back({k});
    if (n_pairs > 0) {
        std::vector<std::tuple<double, Key, Key>> coupled;
        for (Key a : driven)
            for (Key b : driven) {
                if (!(a < b) || !held[a] || !held[b] || a.first > parent_l_max || b.first > parent_l_max) continue;
                double e = INF;
                for (const auto& [d, ea] : best[a])
                    if (best[b].count(d)) e = std::min(e, std::max(ea, best[b].at(d)));
                for (const auto& t : direct)
                    if (std::set<Key>{t.pair[0], t.pair[1]} == std::set<Key>{a, b}) e = std::min(e, energy(t));
                if (std::isfinite(e)) coupled.push_back({e, a, b});
            }
        std::sort(coupled.begin(), coupled.end());
        std::printf("\npairs of parents, each held alone, by the lowest E_th or E_eq of a triplet coupling them:\n");
        for (size_t i=0 ; i<coupled.size() && i<20 ; i++) {
            const auto& [e, a, b] = coupled[i];
            std::printf("  (%d,%+d) (%d,%+d)  %.3e%s\n", a.first, a.second, b.first, b.second, e,
                        int(i) < n_pairs ? "  <- ladder" : "");
            if (int(i) < n_pairs) sets.push_back({a, b});
        }
    }

    for (auto ps_v : sets) {
        std::sort(ps_v.begin(), ps_v.end(), [&](Key a, Key b) {          // slow parent first
            return std::abs(efs.at(a).gamma) < std::abs(efs.at(b).gamma); });
        const std::set<Key> parents(ps_v.begin(), ps_v.end());
        std::string dir = name_of(ps_v[0]);
        if (ps_v.size() == 2) dir += "_" + name_of(ps_v[1]);
        const auto steps = ladder(candidates(parents, only_eeq), parents, n_max);
        if (steps.empty()) {
            std::printf("  %-22s no candidate triplet within the cut, nothing written\n", dir.c_str());
            continue;
        }
        if (only_eeq) {
            const auto plain = ladder(candidates(parents, false), parents, n_max);
            if (plain.size() == steps.size()
                && std::equal(steps.begin(), steps.end(), plain.begin(), [&](const Step& a, const Step& b) {
                       return keys_of(a.chosen) == keys_of(b.chosen); })) {
                std::printf("  %-22s same ladder without --only-eeq, nothing written\n", dir.c_str());
                continue;
            }
        }
        std::filesystem::create_directories(std::filesystem::path(out) / dir);
        for (const Step& s : steps) {
            std::vector<RadialTriplet> net = s.chosen;
            if (closure) {
                ModeMap sub;
                for (Key k : parents) sub.emplace(k, efs.at(k));
                for (Key k : s.daughters) sub.emplace(k, efs.at(k));
                net = enumerate_triplets(sub, closure_cut, l_max);
                const auto have = keys_of(net);
                const std::set<TripletKey> in(have.begin(), have.end());
                for (const auto& t : s.chosen)                // a closure cut below the cut
                    if (!in.count(keys_of({t})[0])) net.push_back(t);
            }
            auto gen = [&](Key k) { return parents.count(k) ? 0 : 1; };
            const auto off = std::remove_if(net.begin(), net.end(), [&](const RadialTriplet& t) {
                return gen(t.pair[0]) != gen(t.pair[1]) || gen(t.sum_mode) == gen(t.pair[0]); });
            const long n_off = long(net.end() - off);
            net.erase(off, net.end());
            stab::kappa_m000(keys_of(net), efs, cache, kap, jobs);
            nonzero(net);

            std::vector<Key> by_id(ps_v.begin(), ps_v.end());
            by_id.insert(by_id.end(), s.daughters.begin(), s.daughters.end());
            char name[32], header[400];
            std::snprintf(name, sizeof name, "d%zu.data", s.daughters.size());
            std::snprintf(header, sizeof header,
                          "# network_build %s  gamma_%s  daughters %zu  cut %g  closure %s  ranking E_th (pairs),"
                          " E_eq (direct)  only_eeq %d  min_gamma_ratio %g  self_coupled %d\n",
                          model.c_str(), gamma_mode.c_str(), s.daughters.size(), cut_dimless,
                          closure ? std::to_string(closure_cut / wdyn).c_str() : "off", int(only_eeq),
                          min_gamma_ratio, n_self);
            write_network(std::filesystem::path(out) / dir / name, header, s, net, by_id, parents, efs, kap);
            std::printf("  %-22s %-9s %3zu selecting, %4zu triplets (%ld off-generation dropped)\n",
                        dir.c_str(), name, s.chosen.size(), net.size(), n_off);
        }
    }
    return 0;
}
