#include "stability.hpp"

#include "csv.hpp"
#include "numeric.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>

namespace stab {

double mu (double kappa, double omega, double delta, double gamma) {
    return std::abs(omega * kappa) / std::hypot(delta, gamma);
}

double parametric_growth_rate (double kappa, double omega_b, double omega_c, double q_a) {
    return 2.0 * std::abs(kappa) * std::sqrt(std::abs(omega_b * omega_c)) * std::abs(q_a);
}

double threshold_energy (double kappa, double omega_b, double omega_c,
                         double gamma_b, double gamma_c, double delta, double E_star) {
    if (gamma_b <= 0.0 || gamma_c <= 0.0 || kappa == 0.0) return 0.0;   // no threshold
    const double g_sum = gamma_b + gamma_c;                             // > 0 by here
    return gamma_b * gamma_c / (4.0 * kappa * kappa * std::abs(omega_b * omega_c))
           * (1.0 + delta * delta / (g_sum * g_sum)) * E_star;
}

double equilibrium_energy (double kappa, double omega_b, double omega_c,
                           double gamma_a, double gamma_b, double gamma_c,
                           double delta, double E_star) {
    if (gamma_b <= 0.0 || gamma_c <= 0.0 || kappa == 0.0) return 0.0;
    const double g_sum = gamma_a + gamma_b + gamma_c;   // the only difference from E_th
    if (g_sum == 0.0) return std::numeric_limits<double>::infinity();
    return gamma_b * gamma_c / (4.0 * kappa * kappa * std::abs(omega_b * omega_c))
           * (1.0 + delta * delta / (g_sum * g_sum)) * E_star;
}

double threshold_energy_ceiling (double kappa, double omega_b, double omega_c,
                                 double delta, double E_star) {
    if (kappa == 0.0) return std::numeric_limits<double>::infinity();
    return delta * delta / (16.0 * kappa * kappa * std::abs(omega_b * omega_c)) * E_star;
}

std::array<double, 3> equilibrium_energies (double kappa, std::array<double, 3> omega,
                                            std::array<double, 3> gamma,
                                            double delta, double E_star) {
    const double nan = std::numeric_limits<double>::quiet_NaN();

    std::array<double, 3> ratio;                 // w_i / gamma_i
    for (int i=0 ; i<3 ; i++) {
        if (gamma[i] == 0.0) return {nan, nan, nan};
        ratio[i] = omega[i] / gamma[i];
        if (!std::isfinite(ratio[i])) return {nan, nan, nan};
    }
    const bool s = ratio[0] > 0.0;
    if ((ratio[1] > 0.0) != s || (ratio[2] > 0.0) != s) return {nan, nan, nan};

    const double E_a = equilibrium_energy(kappa, omega[1], omega[2],
                                          gamma[0], gamma[1], gamma[2], delta, E_star);
    return {E_a, E_a * ratio[1] / ratio[0], E_a * ratio[2] / ratio[0]};
}

eig::Array<bool, eig::Dynamic, 1> nonadiabatic_shell (const Star& s, double omega) {
    return std::abs(omega) * s.t_thermal() < NONADIABATIC_CYCLES;
}

double nonadiabatic_fraction (const KappaResult& res, const Star& s, double omega) {
    const eig::ArrayXd t_th = num::interp(res.r, s.r, s.t_thermal());
    const eig::ArrayXd w = res.dkappa_dr.abs();
    const double total = num::trapz(w, res.r);
    if (total == 0.0) return 0.0;
    const eig::ArrayXd inside =
        (std::abs(omega) * t_th < NONADIABATIC_CYCLES).select(w, 0.0);
    return num::trapz(inside, res.r) / total;
}

const char* channel_name (Channel c) {
    switch (c) {
        case Channel::parametric:  return "parametric";
        case Channel::direct_sum:  return "direct-sum";
        case Channel::direct_diff: return "direct-diff";
        case Channel::all_driven:  return "all-driven";
        case Channel::inactive:    return "inactive";
        case Channel::all_damped:  return "all-damped";
    }
    return "inactive";
}

Channel channel (double gamma_a, double gamma_b, double gamma_c) {
    const bool a = gamma_a < 0.0;                // parents
    const int n_parents = int(a) + int(gamma_b < 0.0) + int(gamma_c < 0.0);
    switch (n_parents) {
        case 3:  return Channel::all_driven;
        case 2:  return a ? Channel::direct_diff : Channel::direct_sum;
        case 1:  return a ? Channel::parametric  : Channel::inactive;
        default: return Channel::all_damped;
    }
}

int daughter_slot (double gamma_a, double gamma_b, double gamma_c) {
    switch (channel(gamma_a, gamma_b, gamma_c)) {
        case Channel::direct_sum:  return 0;                        // a is the daughter
        case Channel::direct_diff: return gamma_b > 0.0 ? 1 : 2;     // whichever is damped
        default:                   return -1;
    }
}

std::vector<TripletObservables> observables (const RadialTriplet& t, const ModeMap& efs,
                                             bool m000) {
    const auto keys = t.keys();
    const Eigenfunction& a = efs.at(keys[0]);
    const Eigenfunction& b = efs.at(keys[1]);
    const Eigenfunction& c = efs.at(keys[2]);
    const double E_star = a.starptr->E_star();
    const std::array<double, 3> gam = {a.gamma, b.gamma, c.gamma};

    const std::vector<std::array<int, 3>> ms =
        m000 ? std::vector<std::array<int, 3>>{{0, 0, 0}} : t.m_combinations();
    const auto [kappas, refine] = kappa_all_m(a, b, c, ms);

    std::vector<TripletObservables> out;
    out.reserve(ms.size());
    for (size_t q=0 ; q<ms.size() ; q++) {
        const double k = kappas[q];
        TripletObservables o;
        o.ms     = ms[q];
        o.kappa  = k;
        o.mu     = {mu(k, t.omega[0], t.delta, gam[0]),
                    mu(k, t.omega[1], t.delta, gam[1]),
                    mu(k, t.omega[2], t.delta, gam[2])};
        o.E_threshold   = threshold_energy(k, t.omega[1], t.omega[2],
                                           gam[1], gam[2], t.delta, E_star);
        o.E_equilibrium = equilibrium_energy(k, t.omega[1], t.omega[2],
                                             gam[0], gam[1], gam[2], t.delta, E_star);
        o.refine = refine;
        out.push_back(o);
    }
    return out;
}

std::vector<Row> build_rows (const std::vector<RadialTriplet>& triplets, const ModeMap& efs,
                             const Star& s, bool m000, int jobs) {
    const double E_star = s.E_star();
    std::vector<std::vector<TripletObservables>> per(triplets.size());

    #pragma omp parallel for schedule(dynamic, 16) num_threads(jobs) if (jobs > 1)
    for (long i=0 ; i<long(triplets.size()) ; i++)
        per[size_t(i)] = observables(triplets[size_t(i)], efs, m000);

    std::vector<Row> rows;
    for (size_t i=0 ; i<triplets.size() ; i++) {
        const RadialTriplet& t = triplets[i];
        const auto l = t.ls(), n = t.ns();
        const auto keys = t.keys();
        const double ga = efs.at(keys[0]).gamma;
        const double gb = efs.at(keys[1]).gamma;
        const double gc = efs.at(keys[2]).gamma;

        for (const TripletObservables& o : per[i]) {
            Row r;
            r.l_a = l[0]; r.n_a = n[0]; r.m_a = o.ms[0];
            r.l_b = l[1]; r.n_b = n[1]; r.m_b = o.ms[1];
            r.l_c = l[2]; r.n_c = n[2]; r.m_c = o.ms[2];
            r.omega_a = t.omega[0]; r.omega_b = t.omega[1]; r.omega_c = t.omega[2];
            r.gamma_a = ga;         r.gamma_b = gb;         r.gamma_c = gc;
            r.delta = t.delta;
            r.kappa = o.kappa;
            r.mu_a = o.mu[0]; r.mu_b = o.mu[1]; r.mu_c = o.mu[2];
            r.mu_max = o.mu_max();
            r.E_th_over_E_star = o.E_threshold / E_star;
            r.E_eq_over_E_star = o.E_equilibrium / E_star;
            r.E_th_ceiling_over_E_star =
                threshold_energy_ceiling(o.kappa, t.omega[1], t.omega[2], t.delta);
            r.detuning_dominated = std::abs(t.delta) > std::abs(ga);
            r.refine = o.refine;
            rows.push_back(r);
        }
    }
    // stable, so equal mu_max keeps enumeration order and diffs stay clean
    std::stable_sort(rows.begin(), rows.end(),
                     [](const Row& x, const Row& y) { return x.mu_max > y.mu_max; });
    return rows;
}

const std::array<const char*, 26> ROW_COLUMNS = {
    "l_a", "n_a", "m_a", "l_b", "n_b", "m_b", "l_c", "n_c", "m_c",
    "omega_a", "omega_b", "omega_c", "gamma_a", "gamma_b", "gamma_c",
    "delta", "kappa", "mu_a", "mu_b", "mu_c", "mu_max",
    "E_th_over_E_star", "E_eq_over_E_star", "E_th_ceiling_over_E_star",
    "detuning_dominated", "refine"
};

void write_rows (const std::filesystem::path& p, const std::vector<Row>& rows, bool force) {
    csv::Writer w(p, {ROW_COLUMNS.begin(), ROW_COLUMNS.end()}, force);
    for (const Row& r : rows)
        w.row(r.l_a, r.n_a, r.m_a, r.l_b, r.n_b, r.m_b, r.l_c, r.n_c, r.m_c,
              r.omega_a, r.omega_b, r.omega_c, r.gamma_a, r.gamma_b, r.gamma_c,
              r.delta, r.kappa, r.mu_a, r.mu_b, r.mu_c, r.mu_max,
              r.E_th_over_E_star, r.E_eq_over_E_star, r.E_th_ceiling_over_E_star,
              r.detuning_dominated, r.refine);
}

Leg leg_of (Key a, Key b, Key c) {
    Leg l{a, b, c};
    std::sort(l.begin(), l.end());
    return l;
}

namespace {

std::string cache_key (const Leg& k) {
    std::string s;
    for (int i = 0; i < 3; ++i)
        s += (i ? ";" : "") + std::to_string(k[i].first) + "," + std::to_string(k[i].second);
    return s;
}

}  // namespace

KappaCache load_kappa_cache (const std::filesystem::path& p) {
    KappaCache out;
    std::ifstream f(p);
    if (!f) return out;
    std::string line;
    std::getline(f, line);                       // header
    while (std::getline(f, line)) {
        const size_t t1 = line.find('\t'), t2 = line.find('\t', t1 + 1);
        if (t1 == std::string::npos || t2 == std::string::npos) continue;
        Leg k{};
        const std::string ks = line.substr(0, t1);
        size_t at = 0;
        for (int i = 0; i < 3; ++i) {
            const size_t comma = ks.find(',', at), semi = ks.find(';', at);
            k[i] = {std::stoi(ks.substr(at, comma - at)), std::stoi(ks.substr(comma + 1))};
            at = semi == std::string::npos ? ks.size() : semi + 1;
        }
        out[k] = {std::stod(line.substr(t1 + 1, t2 - t1 - 1)), std::stoi(line.substr(t2 + 1))};
    }
    return out;
}

void kappa_m000 (const std::vector<Leg>& legs, const ModeMap& efs,
                 const std::filesystem::path& cache, KappaCache& kap, int jobs) {
    std::vector<Leg> todo;
    for (const Leg& l : legs) if (!kap.count(l)) todo.push_back(l);
    std::sort(todo.begin(), todo.end());
    todo.erase(std::unique(todo.begin(), todo.end()), todo.end());
    if (todo.empty()) return;

    if (cache.has_parent_path()) std::filesystem::create_directories(cache.parent_path());
    const bool fresh = !std::filesystem::exists(cache) || std::filesystem::file_size(cache) == 0;
    std::ofstream fh(cache, std::ios::app);
    if (fresh) fh << "key\tkappa\trefine\n";
#pragma omp parallel for schedule(dynamic, 8) num_threads(jobs) if (jobs > 1)
    for (long i = 0; i < long(todo.size()); ++i) {
        const Leg& k = todo[size_t(i)];
        const KappaResult r = kappa_abc(efs.at(k[0]), efs.at(k[1]), efs.at(k[2]), {0, 0, 0});
#pragma omp critical
        {
            kap[k] = {r.kappa, r.refine};
            fh << cache_key(k) << '\t' << csv::fmt(r.kappa) << '\t' << r.refine << '\n';
            fh.flush();
        }
    }
}

}  // namespace stab
