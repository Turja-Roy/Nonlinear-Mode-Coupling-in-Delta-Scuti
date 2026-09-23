#include "amplitude.hpp"

#include <boost/numeric/odeint.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace amp {

namespace odeint = boost::numeric::odeint;

static const std::complex<double> I(0.0, 1.0);

Network::Network (std::vector<Mode> modes, std::vector<Triplet> triplets)
    : modes_(std::move(modes)), triplets_(std::move(triplets)) {
    const int n = size();
    omega_.resize(n);
    gamma_.resize(n);
    for (int i=0 ; i<n ; i++) {
        omega_[i] = modes_[i].omega;
        gamma_[i] = modes_[i].gamma;
    }
    for (const Triplet& tr : triplets_) {
        for (int j : {tr.a, tr.b, tr.c})
            if (j < 0 || j >= n) throw std::invalid_argument("amp::Network: mode index");
        // a == b or a == c would need w_a = w_a + w_other: not a resonance.
        if (tr.a == tr.b || tr.a == tr.c)
            throw std::invalid_argument("amp::Network: sum mode repeated in its own pair");
        delta_.push_back(detuning(tr));
    }
}

double Network::detuning (const Triplet& tr) const {
    return omega_[tr.b] + omega_[tr.c] - omega_[tr.a];
}

/* The three equations of the header, with `phase` carrying exp(-i Delta t) in
   the A equations and 1 in the q equations. The pair equations take its
   conjugate, because they are driven by conj() of a pair partner. */
void Network::add_triplet (const Triplet& tr, const State& q, State& out,
                           std::complex<double> phase) const {
    const int a = tr.a, b = tr.b, c = tr.c;
    const std::complex<double> k = tr.kappa, ph = std::conj(phase);

    if (b == c) {                                // a -> d + d, one pair equation
        out[a] +=       I * omega_[a] * k  * q[b] * q[b] * phase;
        out[b] += 2.0 * I * omega_[b] * k * q[a] * std::conj(q[b]) * ph;
        return;
    }
    out[a] += 2.0 * I * omega_[a] * k * q[b] * q[c] * phase;
    out[b] += 2.0 * I * omega_[b] * k * q[a] * std::conj(q[c]) * ph;
    out[c] += 2.0 * I * omega_[c] * k * q[a] * std::conj(q[b]) * ph;
}

State Network::dq_dt (double t, const State& q) const {
    State out(q.size());
    dq_dt(t, q, out);
    return out;
}

State Network::dA_dt (double t, const State& A) const {
    State out(A.size());
    dA_dt(t, A, out);
    return out;
}

void Network::dq_dt (double t, const State& q, State& out) const {
    (void)t;
    for (int i=0 ; i<size() ; i++) out[i] = -(I * omega_[i] + gamma_[i]) * q[i];
    for (const Triplet& tr : triplets_) add_triplet(tr, q, out, 1.0);
}

void Network::dA_dt (double t, const State& A, State& out) const {
    for (int i=0 ; i<size() ; i++) out[i] = -gamma_[i] * A[i];
    for (size_t j=0 ; j<triplets_.size() ; j++)
        add_triplet(triplets_[j], A, out, std::exp(-I * delta_[j] * t));
}

eig::ArrayXd Network::energy (const State& y) const {
    eig::ArrayXd e(size());
    for (int i=0 ; i<size() ; i++) e[i] = std::norm(y[i]);
    return e;
}

/* One triplet takes f quanta out of the sum mode and puts f into each pair
   mode, so N_a + N_b and N_b - N_c are the invariants; in a network the sums
   run over every triplet a mode belongs to. */
eig::ArrayXd Network::action (const State& y) const {
    return energy(y) / omega_;
}

State Network::to_q (double t, const State& A) const {
    State q(A.size());
    for (int i=0 ; i<size() ; i++) q[i] = A[i] * std::exp(-I * omega_[i] * t);
    return q;
}

namespace {
struct EnergyCap {};                             // thrown to stop the run
}

Solution Network::integrate (State y0, double t_end, Options opt) const {
    if (int(y0.size()) != size())
        throw std::invalid_argument("amp::Network::integrate: y0 has the wrong length");

    if (opt.atol == 0.0) {
        double mx = 0.0;
        for (const auto& z : y0) mx = std::max(mx, std::abs(z));
        opt.atol = 1e-14 * mx;
    }

    auto system = [&](const State& y, State& dydt, double t) {
        if (opt.slow) dA_dt(t, y, dydt);
        else              dq_dt(t, y, dydt);
        for (int i : opt.frozen) dydt[i] = 0.0;
    };

    std::vector<double> times(size_t(opt.n_out));
    for (int i=0 ; i<opt.n_out ; i++) times[size_t(i)] = t_end * i / (opt.n_out - 1);

    Solution sol;
    double e_prev = -1.0;
    auto observe = [&](const State& y, double t) {
        if (opt.e_max > 0.0) {
            double e = 0.0;
            for (const auto& z : y) e += std::norm(z);
            if (e_prev >= 0.0 && e_prev <= opt.e_max && e > opt.e_max) throw EnergyCap{};
            e_prev = e;
        }
        sol.t.push_back(t);
        sol.y.push_back(y);
    };

    auto stepper = odeint::make_controlled<odeint::runge_kutta_fehlberg78<State>>(
        opt.atol, opt.rtol);
    try {
        odeint::integrate_times(stepper, system, y0, times.begin(), times.end(),
                                t_end / (10.0 * opt.n_out), observe);
    } catch (const EnergyCap&) {
    }
    return sol;
}

Network three_mode (std::array<double, 3> omega, std::array<double, 3> gamma, double kappa,
                    std::array<std::string, 3> names) {
    std::vector<Mode> modes;
    for (int i=0 ; i<3 ; i++) modes.push_back(Mode{names[i], omega[i], gamma[i]});
    return Network(std::move(modes), {Triplet{0, 1, 2, kappa}});
}

Network self_coupled (double omega_a, double omega_d, double gamma_a, double gamma_d,
                      double kappa) {
    return Network({Mode{"a", omega_a, gamma_a}, Mode{"d", omega_d, gamma_d}},
                   {Triplet{0, 1, 1, kappa}});
}

Network from_triplets (const std::vector<RadialTriplet>& triplets, const ModeMap& efs,
                       const std::vector<std::array<int, 3>>& ms_list,
                       const std::map<Key, double>& gamma_override) {
    std::map<Key, int> index;
    std::vector<Mode> modes;
    std::vector<Triplet> out;

    auto node = [&](const Key& k) {              // one node per (l, n_pg)
        auto it = index.find(k);
        if (it != index.end()) return it->second;
        const auto ov = gamma_override.find(k);
        const std::string name = "(" + std::to_string(k.first) + "," +
                                 (k.second >= 0 ? "+" : "") + std::to_string(k.second) + ")";
        modes.push_back(Mode{name, efs.at(k).omega,
                             ov != gamma_override.end() ? ov->second : efs.at(k).gamma});
        return index.emplace(k, int(modes.size()) - 1).first->second;
    };

    for (size_t q=0 ; q<triplets.size() ; q++) {
        const RadialTriplet& t = triplets[q];
        const auto keys = t.keys();              // sum mode first
        const auto ms = t.m_combinations();
        const auto [ks, refine] = kappa_all_m(efs.at(keys[0]), efs.at(keys[1]),
                                              efs.at(keys[2]), ms);
        (void)refine;

        size_t pick = 0;
        if (q < ms_list.size()) {
            pick = size_t(std::find(ms.begin(), ms.end(), ms_list[q]) - ms.begin());
            if (pick >= ms.size())
                throw std::invalid_argument("amp::from_triplets: m not in triplet");
        } else {
            for (size_t j=1 ; j<ks.size() ; j++)
                if (std::abs(ks[j]) > std::abs(ks[pick])) pick = j;
        }
        // separate statements: node() interns a new mode, so the order matters
        const int a = node(keys[0]), b = node(keys[1]), c = node(keys[2]);
        out.push_back(Triplet{a, b, c, ks[pick]});
    }
    return Network(std::move(modes), std::move(out));
}

}  // namespace amp
