/* MW25 four-mode mixed systems on the wide net. Two branches, two tables.

       daughter   a + b -> c,  then c -> d + d      resonant with omega_c/2
       parent     a + b -> c,  and a -> d + d,  b -> d + d   at omega_a/2 ~ omega_b/2

   Parent means self-excited (gamma < 0). In the daughter branch c is damped and
   born from a and b, so its own decay products d are granddaughters; nothing
   there absorbs the parents' kappa-mechanism flux, so that branch can only
   answer whether the combination frequency goes parametrically unstable.

   The parent branch is MW25's own configuration: one damped daughter pair, here
   self-coupled, sits at half of both parent frequencies and drains them
   directly. That is what bounds the parents, so it decides whether delta Sct
   parents can be parametrically saturated at all.

   Ranking runs on kappa at m = (0, 0, 0); the top of the ranking then gets the
   full m treatment. */

#include "csv.hpp"
#include "stability.hpp"

#include <CLI/CLI.hpp>

#include <chrono>
#include <cstdio>
#include <limits>
#include <set>

namespace {

constexpr double Q_PARENT = 1e-6;                // MW23 Fig. 6 reference parent amplitude
using stab::Leg;

}  // namespace

int main(int argc, char** argv) {
    CLI::App app{"MW25 four-mode mixed systems on the wide net"};
    std::string model = "models/dsct_M2.0", detail_dir = "detail_wide";
    std::string inlist = "gyre_ad_wide.in", out = "out/four_mode_dsct_M2.0.csv";
    std::vector<std::string> nad = {"summary_nad_wide.h5", "summary_nad_wide_hi.h5"};
    int top_full_m = 200, jobs = 4;
    app.add_option("--model", model);
    app.add_option("--detail-dir", detail_dir);
    app.add_option("--inlist", inlist);
    app.add_option("--nad", nad)->expected(-1);
    app.add_option("--out", out);
    app.add_option("--top-full-m", top_full_m);
    app.add_option("-j,--jobs", jobs);
    CLI11_PARSE(app, argc, argv);

    const auto t0 = std::chrono::steady_clock::now();
    auto secs = [&] {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };

    Model m = load_model(model, detail_dir, inlist, nad);
    ModeMap efs;
    for (auto& [k, e] : m.eigfuncs) if (std::isfinite(e.gamma)) efs.emplace(k, e);
    const double cut = DETUNING_CUT_DIMLESS * m.star->omega_dyn();

    std::vector<Key> keys;
    for (const auto& [k, e] : efs) keys.push_back(k);
    std::sort(keys.begin(), keys.end(),
              [&](Key a, Key b) { return efs.at(a).omega < efs.at(b).omega; });
    const size_t N = keys.size();
    std::vector<double> w(N), g(N);
    std::vector<int> L(N);
    for (size_t i = 0; i < N; ++i) {
        w[i] = efs.at(keys[i]).omega;
        g[i] = efs.at(keys[i]).gamma;
        L[i] = keys[i].first;
    }
    std::vector<int> driven;
    std::map<int, std::vector<int>> by_l;        // per l, sorted by omega like `keys`
    for (size_t i = 0; i < N; ++i) {
        if (g[i] < 0.0) driven.push_back(int(i));
        by_l[L[i]].push_back(int(i));
    }
    std::printf("%zu modes with gamma, %zu driven, cut %.3f c/d\n",
                N, driven.size(), cut / CD);

    struct Tri { int ia, ib, ic, sgn; };
    std::vector<Tri> direct;
    for (size_t i = 0; i < driven.size(); ++i)
        for (size_t jj = i; jj < driven.size(); ++jj) {
            const int ia = driven[i], ib = driven[jj];
            for (auto [wc, sgn] : {std::pair<double, int>{w[ia] + w[ib], +1},
                                   std::pair<double, int>{std::abs(w[ia] - w[ib]), -1}}) {
                if (wc <= 0.05 * CD) continue;
                for (const auto& [lc, idx] : by_l) {
                    if (!satisfies_selection_rules(L[ia], L[ib], lc)) continue;
                    for (int k : idx)
                        if (std::abs(w[k] - wc) < cut && g[k] > 0.0)
                            direct.push_back({ia, ib, k, sgn});
                }
            }
        }

    struct Quad { int ia, ib, ic, id, sgn; };
    std::vector<Quad> quads, par_quads;
    for (const Tri& t : direct) {
        if (L[t.ic] % 2) continue;
        for (const auto& [ld, idx] : by_l) {
            if (2 * ld < L[t.ic]) continue;
            for (int k : idx)
                if (std::abs(2 * w[k] - w[t.ic]) < cut && g[k] > 0.0)
                    quads.push_back({t.ia, t.ib, t.ic, k, t.sgn});
        }
    }
    // MW25's branch: one self-coupled pair at half of *both* parent frequencies.
    for (const Tri& t : direct)
        for (const auto& [ld, idx] : by_l) {
            if (!satisfies_selection_rules(L[t.ia], ld, ld)
                || !satisfies_selection_rules(L[t.ib], ld, ld)) continue;
            for (int k : idx)
                if (std::abs(2 * w[k] - w[t.ia]) < cut
                    && std::abs(2 * w[k] - w[t.ib]) < cut && g[k] > 0.0)
                    par_quads.push_back({t.ia, t.ib, t.ic, k, t.sgn});
        }
    std::printf("%zu direct triples, %zu daughter-branch and %zu parent-branch "
                "candidates [%.0f s]\n", direct.size(), quads.size(), par_quads.size(), secs());
    if (quads.empty() && par_quads.empty()) { std::printf("no candidates\n"); return 1; }

    auto triplet = [&](int sum_slot, int p0, int p1) {
        return RadialTriplet{keys[sum_slot], {keys[p0], keys[p1]},
                             {-w[sum_slot], w[p0], w[p1]}, w[p0] + w[p1] - w[sum_slot]};
    };
    auto leg_of = [&](int i, int j, int k) { return stab::leg_of(keys[i], keys[j], keys[k]); };

    /* The sign assignment does not enter kappa (only omega^2 does), so one
       kappa per unordered (l, n) triple serves every candidate that uses it. */
    std::map<Leg, RadialTriplet> legs;
    for (const Tri& t : direct) {
        int s, p0, p1;
        if (t.sgn > 0)            { s = t.ic; p0 = t.ia; p1 = t.ib; }
        else if (w[t.ia] > w[t.ib]) { s = t.ia; p0 = t.ib; p1 = t.ic; }
        else                        { s = t.ib; p0 = t.ia; p1 = t.ic; }
        legs.emplace(leg_of(t.ia, t.ib, t.ic), triplet(s, p0, p1));
    }
    for (const Quad& q : quads)
        legs.emplace(leg_of(q.ic, q.id, q.id), triplet(q.ic, q.id, q.id));
    for (const Quad& q : par_quads)
        for (int ip : {q.ia, q.ib})
            legs.emplace(leg_of(ip, q.id, q.id), triplet(ip, q.id, q.id));
    std::printf("%zu distinct radial triplets to integrate\n", legs.size());

    const auto t1 = std::chrono::steady_clock::now();
    const std::filesystem::path outp(out);
    std::filesystem::path cache_path = outp;
    cache_path.replace_extension(".kappa_cache.tsv");
    stab::KappaCache kap = stab::load_kappa_cache(cache_path);
    std::vector<Leg> todo;
    for (const auto& [k, t] : legs) if (!kap.count(k)) todo.push_back(k);
    if (!kap.empty())
        std::printf("resuming: %zu triplets cached, %zu to go\n", kap.size(), todo.size());
    stab::kappa_m000(todo, efs, cache_path, kap, jobs);
    std::printf("kappa done, %zu triplets on %d thread(s) [%.0f s]\n", kap.size(), jobs,
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count());

    const double INF = std::numeric_limits<double>::infinity();
    struct DRow { Quad q; double kd, kp; int rd, rp; double dd, dp, q_c, gpar, eth; };
    std::vector<DRow> rows;
    for (const Quad& q : quads) {
        const auto [kd, rd] = kap.at(leg_of(q.ia, q.ib, q.ic));
        const auto [kp, rp] = kap.at(leg_of(q.ic, q.id, q.id));
        const double wa = w[q.ia], wb = w[q.ib], wc = w[q.ic], wd = w[q.id];
        // The sum slot carries the negative sign: the sum branch puts c there,
        // the difference branch the higher-frequency of a, b.
        const double dd = q.sgn > 0 ? wa + wb - wc
                                    : std::min(wa, wb) + wc - std::max(wa, wb);
        const double dp = 2 * wd - wc;
        /* q_c driven off resonance by the two parents, then Gamma of the c -> d,d
           parametric leg at that q_c. s_c = 2 for distinct parents, 1 for a
           self-coupled pair -- the mixed-network integrator measured the 2
           directly against the bare mu. */
        const double s_c = q.ia != q.ib ? 2.0 : 1.0;
        const double q_c = s_c * std::abs(wc * kd) / std::hypot(dd, g[q.ic]) * Q_PARENT * Q_PARENT;
        rows.push_back({q, kd, kp, rd, rp, dd, dp, q_c, 2 * std::abs(kp) * wd * q_c,
                        stab::threshold_energy(kp, wd, wd, g[q.id], g[q.id], dp)});
    }
    struct PRow { Quad q; double kd, ka, kb; int rd, ra, rb; double dd, dpa, dpb,
                  eth_a, eth_b; };
    std::vector<PRow> par_rows;
    for (const Quad& q : par_quads) {
        const auto [kd, rd] = kap.at(leg_of(q.ia, q.ib, q.ic));
        const auto [ka, ra] = kap.at(leg_of(q.ia, q.id, q.id));
        const auto [kb, rb] = kap.at(leg_of(q.ib, q.id, q.id));
        const double wa = w[q.ia], wb = w[q.ib], wc = w[q.ic], wd = w[q.id];
        const double dd = q.sgn > 0 ? wa + wb - wc
                                    : std::min(wa, wb) + wc - std::max(wa, wb);
        const double dpa = 2 * wd - wa, dpb = 2 * wd - wb;
        /* Not E_c/E_th: what decides this branch is whether a parent at the
           observed amplitude crosses its own parametric threshold, so the
           relevant threshold is the smaller of the two parents'. */
        par_rows.push_back({q, kd, ka, kb, rd, ra, rb, dd, dpa, dpb,
                            stab::threshold_energy(ka, wd, wd, g[q.id], g[q.id], dpa),
                            stab::threshold_energy(kb, wd, wd, g[q.id], g[q.id], dpb)});
    }

    auto ratio = [&](double num, double den) { return den != 0.0 ? num / den : INF; };
    if (!par_rows.empty()) {
        std::stable_sort(par_rows.begin(), par_rows.end(), [&](const PRow& a, const PRow& b) {
            return ratio(Q_PARENT * Q_PARENT, std::min(a.eth_a, a.eth_b))
                 > ratio(Q_PARENT * Q_PARENT, std::min(b.eth_a, b.eth_b));
        });
        std::filesystem::path pout = outp;
        pout.replace_filename(outp.stem().string() + "_parent.csv");
        csv::Writer wcsv(pout, {"branch", "l_a", "n_a", "l_b", "n_b", "l_c", "n_c", "l_d", "n_d",
                              "comb", "f_a", "f_b", "f_c", "f_d",
                              "gamma_a", "gamma_b", "gamma_c", "gamma_d",
                              "delta_direct", "delta_param_a", "delta_param_b",
                              "kappa_direct", "kappa_param_a", "kappa_param_b",
                              "refine_direct", "refine_param_a", "refine_param_b",
                              "E_th_par_a", "E_th_par_b", "E_th_par", "E_par_over_E_th"});
        long above = 0;
        for (const PRow& r : par_rows) {
            const double eth = std::min(r.eth_a, r.eth_b);
            const double rel = ratio(Q_PARENT * Q_PARENT, eth);
            above += rel > 1.0;
            wcsv.row("parent", keys[r.q.ia].first, keys[r.q.ia].second,
                     keys[r.q.ib].first, keys[r.q.ib].second,
                     keys[r.q.ic].first, keys[r.q.ic].second,
                     keys[r.q.id].first, keys[r.q.id].second,
                     r.q.sgn > 0 ? "sum" : "diff",
                     w[r.q.ia] / CD, w[r.q.ib] / CD, w[r.q.ic] / CD, w[r.q.id] / CD,
                     g[r.q.ia], g[r.q.ib], g[r.q.ic], g[r.q.id],
                     r.dd, r.dpa, r.dpb, r.kd, r.ka, r.kb, r.rd, r.ra, r.rb,
                     r.eth_a, r.eth_b, eth, rel);
        }
        std::printf("\n%zu parent-branch rows -> %s\n", par_rows.size(), pout.c_str());
        std::printf("%ld parents above their own parametric threshold at q_parent = %g\n",
                    above, Q_PARENT);
    }

    if (rows.empty()) { std::printf("\ntotal %.0f s\n", secs()); return 0; }
    std::stable_sort(rows.begin(), rows.end(), [&](const DRow& a, const DRow& b) {
        return ratio(a.q_c * a.q_c, a.eth) > ratio(b.q_c * b.q_c, b.eth);
    });
    {
        csv::Writer wcsv(outp, {"branch", "l_a", "n_a", "l_b", "n_b", "l_c", "n_c", "l_d", "n_d",
                              "comb", "f_a", "f_b", "f_c", "f_d",
                              "gamma_a", "gamma_b", "gamma_c", "gamma_d",
                              "delta_direct", "delta_param", "kappa_direct", "kappa_param",
                              "refine_direct", "refine_param", "q_c", "gamma_par",
                              "growth_over_damping", "E_c_over_E_star",
                              "E_th_over_E_star", "E_c_over_E_th"});
        long above = 0, above_bare = 0;
        for (const DRow& r : rows) {
            const double rel = ratio(r.q_c * r.q_c, r.eth);
            above += rel > 1.0;
            above_bare += r.gpar / g[r.q.id] > 1.0;
            wcsv.row("daughter", keys[r.q.ia].first, keys[r.q.ia].second,
                     keys[r.q.ib].first, keys[r.q.ib].second,
                     keys[r.q.ic].first, keys[r.q.ic].second,
                     keys[r.q.id].first, keys[r.q.id].second,
                     r.q.sgn > 0 ? "sum" : "diff",
                     w[r.q.ia] / CD, w[r.q.ib] / CD, w[r.q.ic] / CD, w[r.q.id] / CD,
                     g[r.q.ia], g[r.q.ib], g[r.q.ic], g[r.q.id],
                     r.dd, r.dp, r.kd, r.kp, r.rd, r.rp, r.q_c, r.gpar,
                     r.gpar / g[r.q.id], r.q_c * r.q_c, r.eth, rel);
        }
        std::printf("\n%zu rows -> %s\n", rows.size(), outp.c_str());
        /* Gamma > gamma_d is not the criterion: the parametric threshold carries
           the detuning factor, and |Delta_p|/2gamma_d has median ~89 here, so
           ranking on Gamma/gamma_d overstates the instability by that squared. */
        std::printf("\n%ld candidates above threshold (E_c > E_th) at q_parent = %g; "
                    "%ld if the detuning factor is dropped\n", above, Q_PARENT, above_bare);
    }

    std::printf("\nfull m for the top %d:\n", top_full_m);
    std::filesystem::path fmout = outp;
    fmout.replace_filename(outp.stem().string() + "_fullm.csv");
    csv::Writer fw(fmout, {"l_c", "n_c", "l_d", "n_d", "kappa_m0", "kappa_max", "m_max",
                           "n_m", "refine"});
    std::set<std::pair<Key, Key>> seen;
    int done = 0;
    for (const DRow& r : rows) {
        const std::pair<Key, Key> key{keys[r.q.ic], keys[r.q.id]};
        if (!seen.insert(key).second) continue;
        const RadialTriplet& t = legs.at(leg_of(r.q.ic, r.q.id, r.q.id));
        const auto tk = t.keys();
        const auto ms = t.m_combinations();
        const auto [ks, refine] = kappa_all_m(efs.at(tk[0]), efs.at(tk[1]), efs.at(tk[2]), ms);
        size_t best = 0;
        for (size_t j = 1; j < ks.size(); ++j)
            if (std::abs(ks[j]) > std::abs(ks[best])) best = j;
        fw.row(key.first.first, key.first.second, key.second.first, key.second.second,
               r.kp, ks[best],
               "(" + std::to_string(ms[best][0]) + ", " + std::to_string(ms[best][1])
                   + ", " + std::to_string(ms[best][2]) + ")",
               long(ks.size()), refine);
        if (++done >= top_full_m) break;
    }
    std::printf("-> %s\ntotal %.0f s\n", fmout.c_str(), secs());
    return 0;
}
