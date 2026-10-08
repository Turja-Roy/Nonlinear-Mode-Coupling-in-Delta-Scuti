#include "amplitude.hpp"

#include <arkode/arkode_arkstep.h>
#include <cvode/cvode.h>
#include <nvector/nvector_openmp.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace amp {

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
        // Only b == c (a -> b + b) is a real self-coupling; a == b or a == c
        // would mean a couples to itself and one other mode, not a triplet.
        if (tr.a == tr.b || tr.a == tr.c)
            throw std::invalid_argument("amp::Network: mode repeated across triplet roles");
        delta_.push_back(detuning(tr));
    }
}

double Network::detuning (const Triplet& tr) const {
    return omega_[tr.a] + omega_[tr.b] + omega_[tr.c];
}

/* -------------------------
    The Amplitude Equations
   ------------------------- */
/* The equations of the header, with `phase` carrying exp(i Delta t) in the A
   equations and 1 in the q equations -- the same value for every slot, so
   unlike a sum-slot scheme there is nothing here to get backwards. */
void Network::add_triplet (const Triplet& tr, const State& q, State& out,
                           std::complex<double> phase) const {
    const int a = tr.a, b = tr.b, c = tr.c;
    const std::complex<double> k = tr.kappa;

    if (b == c) {                                // a -> b + b, one pair equation
        out[a] +=       I * omega_[a] * k * std::conj(q[b]) * std::conj(q[b]) * phase;
        out[b] += 2.0 * I * omega_[b] * k * std::conj(q[a]) * std::conj(q[b]) * phase;
        return;
    }
    out[a] += 2.0 * I * omega_[a] * k * std::conj(q[b]) * std::conj(q[c]) * phase;
    out[b] += 2.0 * I * omega_[b] * k * std::conj(q[a]) * std::conj(q[c]) * phase;
    out[c] += 2.0 * I * omega_[c] * k * std::conj(q[a]) * std::conj(q[b]) * phase;
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
    for (int i=0 ; i<size() ; i++) out[i] = -(I * omega_[i] + gamma_[i]) * q[i];
    add_triplets(t, q, out, false);
}

void Network::dA_dt (double t, const State& A, State& out) const {
    for (int i=0 ; i<size() ; i++) out[i] = -gamma_[i] * A[i];
    add_triplets(t, A, out, true);
}

/* Triplets are split over OpenMP threads once there are enough of them to pay
   for the threads. Two triplets sharing a mode both write out[a], so each
   thread sums into its own buffer and the buffers are added at the end.
   Called from inside an outer parallel loop (stability.cpp, four_mode_search)
   this runs on one thread, since OpenMP nesting is off by default. */
void Network::add_triplets (double t, const State& y, State& out, bool slow) const {
    const size_t nt = triplets_.size();
    auto phase = [&](size_t j) {
        return slow ? std::exp(I * delta_[j] * t) : std::complex<double>(1.0);
    };
    if (nt <= 64) {                              // small network: threads cost more than they save
        for (size_t j=0 ; j<nt ; j++) add_triplet(triplets_[j], y, out, phase(j));
        return;
    }
    // ponytail: one buffer per thread per call; per-mode ownership if this shows in profiles
    #pragma omp parallel
    {
        State local(out.size(), 0.0);
        #pragma omp for nowait
        for (size_t j=0 ; j<nt ; j++) add_triplet(triplets_[j], y, local, phase(j));
        #pragma omp critical
        for (size_t i=0 ; i<out.size() ; i++) out[i] += local[i];
    }
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

/* --------------------
    The ODE Integrator
   -------------------- */
namespace {

struct SysParams {
    const Network* net;
    const Options* opt;
    State y, dydt;                               // scratch, so the rhs never allocates
};

/* The solvers work with raw doubles, y[2n] = {Re q_0, Im q_0, Re q_1, ...};
   std::complex<double> has exactly that layout, so a cast converts. */
std::complex<double>* as_complex (double* y) {
    return reinterpret_cast<std::complex<double>*>(y);
}

// Threads for SUNDIALS' vector operations: they are a few flops per element,
// so threading only pays once the state is large.
int threads (int n) {
#ifdef _OPENMP
    // ponytail: fixed cut; measure on the cluster once networks reach it
    return n > 5000 ? omp_get_max_threads() : 1;
#else
    (void)n;
    return 1;
#endif
}

// Full right-hand side, dA/dt (or dq/dt), with frozen modes held.
void full_rhs (double t, const double* y, double* dydt, SysParams* p) {
    const size_t n = p->y.size();
    std::copy_n(reinterpret_cast<const std::complex<double>*>(y), n, p->y.begin());

    if (p->opt->slow) p->net->dA_dt(t, p->y, p->dydt);
    else              p->net->dq_dt(t, p->y, p->dydt);
    for (int i : p->opt->frozen) p->dydt[size_t(i)] = 0.0;

    std::copy_n(p->dydt.begin(), n, as_complex(dydt));
}

// The linear term alone: -g y for A, -(i w + g) y for q.
std::complex<double> linear_rate (const Mode& m, bool slow) {
    return slow ? std::complex<double>(-m.gamma) : -(I * m.omega + m.gamma);
}

// CVODE: the whole right-hand side.
[[maybe_unused]] int rhs (sunrealtype t, N_Vector y, N_Vector ydot, void* params) {
    full_rhs(t, N_VGetArrayPointer(y), N_VGetArrayPointer(ydot), static_cast<SysParams*>(params));
    return 0;
}

// ARKODE IMEX, implicit part: the linear (damping) term, where the stiffness is.
[[maybe_unused]] int rhs_implicit (sunrealtype t, N_Vector y, N_Vector ydot, void* params) {
    (void)t;
    const auto* p = static_cast<SysParams*>(params);
    const std::complex<double>* in = as_complex(N_VGetArrayPointer(y));
    std::complex<double>* out = as_complex(N_VGetArrayPointer(ydot));
    for (int i=0 ; i<p->net->size() ; i++) out[i] = linear_rate(p->net->modes()[size_t(i)], p->opt->slow) * in[i];
    for (int i : p->opt->frozen) out[i] = 0.0;
    return 0;
}

// ARKODE IMEX, explicit part: the triplet coupling, i.e. the full rhs minus the linear term.
[[maybe_unused]] int rhs_explicit (sunrealtype t, N_Vector y, N_Vector ydot, void* params) {
    auto* p = static_cast<SysParams*>(params);
    full_rhs(t, N_VGetArrayPointer(y), N_VGetArrayPointer(ydot), p);
    const std::complex<double>* in = as_complex(N_VGetArrayPointer(y));
    std::complex<double>* out = as_complex(N_VGetArrayPointer(ydot));
    for (int i=0 ; i<p->net->size() ; i++) out[i] -= linear_rate(p->net->modes()[size_t(i)], p->opt->slow) * in[i];
    for (int i : p->opt->frozen) out[i] = 0.0;
    return 0;
}

// GSL: the whole right-hand side, GSL's signature.
[[maybe_unused]] int rhs_gsl (double t, const double y[], double dydt[], void* params) {
    full_rhs(t, y, dydt, static_cast<SysParams*>(params));
    return GSL_SUCCESS;
}

}  // namespace

Solution Network::integrate (State y0, double t_end, Options opt) const {
    const int n = size();
    if (int(y0.size()) != n)
        throw std::invalid_argument("amp::Network::integrate: y0 has the wrong length");

    if (opt.atol == 0.0) {
        double mx = 0.0;
        for (const auto& z : y0) mx = std::max(mx, std::abs(z));
        opt.atol = 1e-14 * mx;
    }

    SysParams params{this, &opt, State(size_t(n)), State(size_t(n))};

    // The state vector: 2n doubles (real + imaginary parts).
    SUNContext ctx;
    SUNContext_Create(SUN_COMM_NULL, &ctx);
    N_Vector y = N_VNew_OpenMP(2 * n, threads(n), ctx);
    double* yd = N_VGetArrayPointer(y);
    std::copy_n(y0.begin(), n, as_complex(yd));


    // (1) CVODE, BDF + GMRES: implicit, variable order 1-5, safe when the
    //     damping rates span orders of magnitude (stiff). GMRES solves the
    //     Newton systems matrix-free, so no Jacobian is needed.
    // void* mem = CVodeCreate(CV_BDF, ctx);
    // CVodeInit(mem, rhs, 0.0, y);
    // CVodeSStolerances(mem, opt.rtol, opt.atol);
    // CVodeSetUserData(mem, &params);
    // CVodeSetMaxNumSteps(mem, -1);                // no cap on steps per output interval
    // SUNLinearSolver LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, ctx);
    // CVodeSetLinearSolver(mem, LS, nullptr);

    // (2) CVODE, Adams + fixed point: multistep, variable order 1-12, no
    //     linear solves. Non-stiff only; cheapest when the rates are comparable.
    // void* mem = CVodeCreate(CV_ADAMS, ctx);
    // CVodeInit(mem, rhs, 0.0, y);
    // CVodeSStolerances(mem, opt.rtol, opt.atol);
    // CVodeSetUserData(mem, &params);
    // CVodeSetMaxNumSteps(mem, -1);
    // SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(y, 0, ctx);
    // CVodeSetNonlinearSolver(mem, NLS);

    // (3) ARKODE, IMEX additive Runge-Kutta: the damping term (rhs_implicit)
    //     is treated implicitly, the triplet coupling (rhs_explicit) explicitly.
    // void* mem = ARKStepCreate(rhs_explicit, rhs_implicit, 0.0, y, ctx);
    // ARKodeSStolerances(mem, opt.rtol, opt.atol);
    // ARKodeSetUserData(mem, &params);
    // ARKodeSetMaxNumSteps(mem, -1);
    // ARKodeSetLinear(mem, 0);                  // rhs_implicit is linear in y
    // SUNLinearSolver LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, ctx);
    // ARKodeSetLinearSolver(mem, LS, nullptr);

    // (4) GSL, Prince-Dormand 8(9): explicit adaptive Runge-Kutta, non-stiff.
    //     Initial step t_end / (10 n_out).
    gsl_odeiv2_system sys{rhs_gsl, nullptr, size_t(2 * n), &params};
    /* First trial step well inside the fastest rate. Too large a step overflows
       to NaN, every error comparison against NaN is false, and the driver then
       accepts it and runs on to t_end on garbage. */
    double rate = 0.0;
    for (int i=0 ; i<n ; i++) rate = std::max(rate, std::abs(gamma_[i]));
    for (double d : delta_) rate = std::max(rate, std::abs(d));
    const double h0 = rate > 0.0 ? std::min(t_end / (10.0 * opt.n_out), 0.01 / rate)
                                 : t_end / (10.0 * opt.n_out);
    gsl_odeiv2_driver* drv = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rk8pd, h0, opt.atol, opt.rtol);

    /* The integration loop: n_out outputs evenly spaced over [0, t_end]. With
       a cap, energy is also checked every e-folding of the fastest driven mode,
       so a runaway stops near the crossing instead of overflowing inside one
       long output interval; the crossing state is then the last output. */
    double grow = 0.0;
    for (int i=0 ; i<n ; i++) grow = std::max(grow, -2.0 * gamma_[i]);
    const double check = opt.e_max > 0.0 && grow > 0.0 ? 1.0 / grow : t_end;
    Solution sol;
    State out(static_cast<size_t>(n));
    double t = 0.0, e_prev = -1.0;
    auto energy_now = [&] {
        std::copy_n(as_complex(yd), n, out.begin());
        double e = 0.0;
        for (const auto& z : out) e += std::norm(z);
        return e;
    };
    auto crossed = [&](double e) {
        const bool c = opt.e_max > 0.0 && e_prev >= 0.0 && e_prev <= opt.e_max && e > opt.e_max;
        if (opt.e_max > 0.0) e_prev = e;
        return c;
    };
    bool done = false;
    for (int step=0 ; step<opt.n_out && !done ; step++) {
        const double t_target = t_end * step / (opt.n_out - 1);

        // advance from t to t_target, updating t and y in place; flag < 0 is a failure
        int flag = 0;
        double e = energy_now();
        while (t < t_target && flag == 0) {
            const double t_next = std::min(t_target, t + check);
            // flag = CVode(mem, t_next, y, &t, CV_NORMAL);                     // (1), (2)
            // flag = ARKodeEvolve(mem, t_next, y, &t, ARK_NORMAL);          // (3)
            flag = gsl_odeiv2_driver_apply(drv, &t, t_next, yd) == GSL_SUCCESS ? 0 : -1;  // (4)
            e = energy_now();
            if (!std::isfinite(e)) { flag = -1; break; }
            if (t < t_target && crossed(e)) { done = true; break; }
        }
        if (flag < 0) break;
        sol.t.push_back(t);
        sol.y.push_back(out);
        if (!done && crossed(e)) done = true;
    }

    // SUNLinSolFree(LS);                           // (1), (3)
    // SUNNonlinSolFree(NLS);                    // (2)
    // CVodeFree(&mem);                             // (1), (2)
    // ARKodeFree(&mem);                         // (3)
    gsl_odeiv2_driver_free(drv);              // (4)
    N_VDestroy(y);
    SUNContext_Free(&ctx);
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

}  // namespace amp
