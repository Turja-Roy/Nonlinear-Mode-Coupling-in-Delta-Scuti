#include "triplets.hpp"

#include <algorithm>
#include <cmath>

std::array<int, 3> RadialTriplet::ls() const {
    const auto k = keys();
    return {k[0].first, k[1].first, k[2].first};
}

std::array<int, 3> RadialTriplet::ns() const {
    const auto k = keys();
    return {k[0].second, k[1].second, k[2].second};
}

std::vector<std::array<int, 3>> RadialTriplet::m_combinations() const {
    const auto l = ls();
    return m_combos(l[0], l[1], l[2]);
}

std::array<Mode, 3> RadialTriplet::modes(const ModeMap& efs, std::array<int, 3> ms) const {
    const auto k = keys();
    const auto s = signs();
    std::array<Mode, 3> out;
    for (int i = 0; i < 3; ++i)
        out[i] = Mode{k[i].second, k[i].first, ms[i], s[i], &efs.at(k[i])};
    return out;
}

namespace {

// One mode list per l, sorted by omega. Parallel arrays so the frequency
// column can be binary-searched directly.
struct ByL {
    std::vector<int> n;
    std::vector<double> w;
};

// Distinct (sum-slot l, l, l); the pair is unordered.
std::vector<std::array<int, 3>> sum_slot_assignments(const std::array<int, 3>& lm) {
    std::set<std::pair<int, std::array<int, 2>>> seen;
    std::vector<std::array<int, 3>> out;
    for (int i = 0; i < 3; ++i) {
        std::array<int, 2> rest{lm[(i + 1) % 3], lm[(i + 2) % 3]};
        if (rest[0] > rest[1]) std::swap(rest[0], rest[1]);
        if (seen.insert({lm[i], rest}).second) out.push_back({lm[i], rest[0], rest[1]});
    }
    return out;
}

void match(std::vector<RadialTriplet>& out, const ByL& p, const ByL& q, const ByL& r,
           int l_p, int l_q, int l_r, double cut) {
    for (size_t iq = 0; iq < q.n.size(); ++iq)
        // The pair is unordered, so restrict to the upper triangle when its
        // two members share an l. iq == ir is kept: the self-coupled case.
        for (size_t ir = (l_q == l_r ? iq : 0); ir < r.n.size(); ++ir) {
            const double s = q.w[iq] + r.w[ir];
            const auto lo = std::lower_bound(p.w.begin(), p.w.end(), s - cut);
            const auto hi = std::upper_bound(p.w.begin(), p.w.end(), s + cut);
            for (auto it = lo; it != hi; ++it) {
                const size_t ip = it - p.w.begin();
                out.push_back(RadialTriplet{
                    {l_p, p.n[ip]},
                    {Key{l_q, q.n[iq]}, Key{l_r, r.n[ir]}},
                    {-p.w[ip], q.w[iq], r.w[ir]},
                    s - p.w[ip],
                });
            }
        }
}

}  // namespace

std::vector<RadialTriplet> enumerate_triplets(const ModeMap& efs, double cut, int l_max,
                                              const std::set<Key>* sum_keys) {
    std::map<int, ByL> by_l, sum_l;
    for (int l = 0; l <= l_max; ++l) {
        std::vector<Key> keys;
        for (const auto& [k, ef] : efs) if (k.first == l) keys.push_back(k);
        if (keys.empty()) continue;
        std::sort(keys.begin(), keys.end(),
                  [&](Key a, Key b) { return efs.at(a).omega < efs.at(b).omega; });

        ByL all, sub;
        for (Key k : keys) {
            all.n.push_back(k.second);
            all.w.push_back(efs.at(k).omega);
            if (!sum_keys || sum_keys->count(k)) {
                sub.n.push_back(k.second);
                sub.w.push_back(efs.at(k).omega);
            }
        }
        by_l[l] = std::move(all);
        if (!sub.n.empty()) sum_l[l] = std::move(sub);
    }

    std::vector<RadialTriplet> out;
    for (const auto& lm : l_multisets(l_max))
        for (const auto& [l_p, l_q, l_r] : sum_slot_assignments(lm)) {
            if (!sum_l.count(l_p) || !by_l.count(l_q) || !by_l.count(l_r)) continue;
            match(out, sum_l[l_p], by_l[l_q], by_l[l_r], l_p, l_q, l_r, cut);
        }
    return out;
}

long count_with_m(const std::vector<RadialTriplet>& triplets) {
    long n = 0;
    for (const auto& t : triplets) n += long(t.m_combinations().size());
    return n;
}

TripletSummary summarise(const std::vector<RadialTriplet>& triplets) {
    TripletSummary s{};
    s.n_radial = long(triplets.size());
    if (triplets.empty()) return s;
    s.n_with_m = count_with_m(triplets);

    std::vector<double> frac;
    frac.reserve(triplets.size());
    for (const auto& t : triplets) {
        frac.push_back(std::abs(t.delta) / -t.omega[0]);
        auto l = t.ls();
        std::sort(l.begin(), l.end());
        s.l_combinations[l] += 1;
    }
    std::sort(frac.begin(), frac.end());
    s.min_abs_delta_over_omega = frac.front();
    const size_t h = frac.size() / 2;
    s.median_abs_delta_over_omega =
        frac.size() % 2 ? frac[h] : 0.5 * (frac[h - 1] + frac[h]);
    return s;
}

std::vector<TripletRow> triplet_rows(const std::vector<RadialTriplet>& triplets,
                                     bool with_m) {
    std::vector<TripletRow> out;
    out.reserve(triplets.size());
    for (const auto& t : triplets) {
        const auto l = t.ls(), n = t.ns();
        out.push_back(TripletRow{
            l[0], n[0], l[1], n[1], l[2], n[2],
            t.omega[0], t.omega[1], t.omega[2], t.delta,
            std::abs(t.delta / t.omega[0]),
            with_m ? long(t.m_combinations().size()) : -1,
        });
    }
    return out;
}
