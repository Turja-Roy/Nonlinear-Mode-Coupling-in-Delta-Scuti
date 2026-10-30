#include "amplitude.hpp"

#include <cvode/cvode.h>
#include <nvector/nvector_openmp.h>
#include <nvector/nvector_serial.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#if 0   // the solvers set aside, see Options
#include <cvode/cvode_bandpre.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <chrono>
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

void Network::jacobian (double t, const State& y, bool slow, double* J) const {
    const size_t m = 2 * size_t(size());
    std::fill(J, J + m * m, 0.0);
    auto at = [&](int r, int c) -> double& { return J[size_t(c) * m + size_t(r)]; };
    auto lin = [&](int x, std::complex<double> L) {          // d(L y_x)/dy_x
        at(2*x, 2*x) += L.real();   at(2*x, 2*x+1) -= L.imag();
        at(2*x+1, 2*x) += L.imag(); at(2*x+1, 2*x+1) += L.real();
    };
    auto anti = [&](int x, int u, std::complex<double> W) {  // d(W conj(y_u))/dy_u, into row x
        at(2*x, 2*u) += W.real();   at(2*x, 2*u+1) += W.imag();
        at(2*x+1, 2*u) += W.imag(); at(2*x+1, 2*u+1) -= W.real();
    };
    for (int i=0 ; i<size() ; i++)
        lin(i, slow ? std::complex<double>(-gamma_[i]) : -(I * omega_[i] + gamma_[i]));
    for (size_t j=0 ; j<triplets_.size() ; j++) {
        const Triplet& tr = triplets_[j];
        const std::complex<double> k = tr.kappa * (slow ? std::exp(I * delta_[j] * t) : 1.0);
        auto term = [&](int x, double c, int u, int v) {     // out[x] += c i w_x k conj(y_u) conj(y_v)
            const std::complex<double> C = c * I * omega_[x] * k;
            anti(x, u, C * std::conj(y[size_t(v)]));
            anti(x, v, C * std::conj(y[size_t(u)]));
        };
        if (tr.b == tr.c) {                      // same terms as add_triplet
            term(tr.a, 1.0, tr.b, tr.b);
            term(tr.b, 2.0, tr.a, tr.b);
            continue;
        }
        term(tr.a, 2.0, tr.b, tr.c);
        term(tr.b, 2.0, tr.a, tr.c);
        term(tr.c, 2.0, tr.a, tr.b);
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
    long rhs_calls = 0;
};

/* The solvers work with raw doubles, y[2n] = {Re q_0, Im q_0, Re q_1, ...};
   std::complex<double> has exactly that layout, so a cast converts. */
std::complex<double>* as_complex (double* y) {
    return reinterpret_cast<std::complex<double>*>(y);
}

// Threads for SUNDIALS' vector operations: they are a few flops per element,
// so threading only pays once the state is large. Below that the serial vector,
// not OpenMP with one thread: every OpenMP vector operation still opens a
// parallel region, ~1 us each, dozens per CVODE step -- 10x the whole step on
// a small network.
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
    ++p->rhs_calls;
    const size_t n = p->y.size();
    std::copy_n(reinterpret_cast<const std::complex<double>*>(y), n, p->y.begin());

    if (p->opt->slow) p->net->dA_dt(t, p->y, p->dydt);
    else              p->net->dq_dt(t, p->y, p->dydt);
    for (int i : p->opt->frozen) p->dydt[size_t(i)] = 0.0;

    std::copy_n(p->dydt.begin(), n, as_complex(dydt));
}

// CVODE: the whole right-hand side.
int rhs (sunrealtype t, N_Vector y, N_Vector ydot, void* params) {
    full_rhs(t, N_VGetArrayPointer(y), N_VGetArrayPointer(ydot), static_cast<SysParams*>(params));
    return 0;
}

/* Error weight per mode, not per real component: |A_i| sets the tolerance on
   both Re A_i and Im A_i. Per component, a phase that puts one of them at 0
   asks for atol absolute there, and a real start with Delta != 0 fails its
   first step. */
int ewt (N_Vector y, N_Vector w, void* params) {
    const auto* p = static_cast<SysParams*>(params);
    const double* yd = N_VGetArrayPointer(y);
    double* wd = N_VGetArrayPointer(w);
    for (size_t i = 0; i < p->y.size(); ++i)
        wd[2 * i] = wd[2 * i + 1] = 1.0 / (p->opt->rtol * std::hypot(yd[2 * i], yd[2 * i + 1]) + p->opt->atol);
    return 0;
}

#if 0
// CVODE BDF: the analytic Jacobian, frozen modes' rows zeroed like their rhs.
int jac (sunrealtype t, N_Vector y, N_Vector, SUNMatrix J, void* params, N_Vector, N_Vector, N_Vector) {
    auto* p = static_cast<SysParams*>(params);
    const size_t n = p->y.size(), m = 2 * n;
    std::copy_n(reinterpret_cast<const std::complex<double>*>(N_VGetArrayPointer(y)), n, p->y.begin());
    double* d = SUNDenseMatrix_Data(J);
    p->net->jacobian(t, p->y, p->opt->slow, d);
    for (int i : p->opt->frozen)
        for (size_t c = 0; c < m; ++c) d[c * m + 2 * size_t(i)] = d[c * m + 2 * size_t(i) + 1] = 0.0;
    return 0;
}

// GSL: the whole right-hand side, GSL's signature.
int rhs_gsl (double t, const double y[], double dydt[], void* params) {
    full_rhs(t, y, dydt, static_cast<SysParams*>(params));
    return GSL_SUCCESS;
}
#endif

}  // namespace

Solution Network::integrate (State y0, double t_end, Options opt) const {
    const int n = size();
    if (int(y0.size()) != n)
        throw std::invalid_argument("amp::Network::integrate: y0 has the wrong length");

    SysParams params{this, &opt, State(size_t(n)), State(size_t(n))};

    // The state vector: 2n doubles (real + imaginary parts).
    SUNContext ctx;
    SUNContext_Create(SUN_COMM_NULL, &ctx);
    N_Vector y = threads(n) > 1 ? N_VNew_OpenMP(2 * n, threads(n), ctx) : N_VNew_Serial(2 * n, ctx);
    double* yd = N_VGetArrayPointer(y);
    std::copy_n(y0.begin(), n, as_complex(yd));

    /* First trial step well inside the fastest rate. Too large a step overflows
       to NaN, every error comparison against NaN is false, and the driver then
       accepts it and runs on to t_end on garbage. */
    double rate = 0.0;
    for (int i=0 ; i<n ; i++) rate = std::max(rate, std::abs(gamma_[i]));
    for (double d : delta_) rate = std::max(rate, std::abs(d));
    const double h0 = rate > 0.0 ? std::min(t_end / (10.0 * opt.n_out), 0.01 / rate)
                                 : t_end / (10.0 * opt.n_out);

    void* mem = CVodeCreate(CV_ADAMS, ctx);
    CVodeInit(mem, rhs, 0.0, y);
    CVodeSetUserData(mem, &params);
    CVodeWFtolerances(mem, ewt);
    CVodeSetInitStep(mem, h0);
    CVodeSetMaxNumSteps(mem, -1);                // no cap on steps per output interval
    SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(y, 0, ctx);
    CVodeSetNonlinearSolver(mem, NLS);
#if 0   // instead of the Adams lines above; GSL also restarts every chunk, see the loop
    gsl_odeiv2_system sys{rhs_gsl, nullptr, size_t(2 * n), &params};
    gsl_odeiv2_driver* drv = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, h0, opt.atol, opt.rtol);
    // BDF, dense
    mem = CVodeCreate(CV_BDF, ctx);
    SUNMatrix J = SUNDenseMatrix(2 * n, 2 * n, ctx);
    SUNLinearSolver LS = SUNLinSol_Dense(y, J, ctx);
    CVodeSetLinearSolver(mem, LS, J);
    CVodeSetJacFn(mem, jac);
    CVodeSetStabLimDet(mem, SUNTRUE);            // drops order when BDF 3-5 meet an oscillation
    // BDF, Krylov: SPGMR, CVBANDPRE preconditioner of half-bandwidth `band`
    LS = SUNLinSol_SPGMR(y, SUN_PREC_LEFT, 0, ctx);
    CVodeSetLinearSolver(mem, LS, nullptr);
    CVBandPrecInit(mem, 2 * n, band, band);
#endif

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
    const auto clock0 = std::chrono::steady_clock::now();
    auto out_of_time = [&] {
        return opt.max_seconds > 0.0 && std::chrono::duration<double>(
                   std::chrono::steady_clock::now() - clock0).count() > opt.max_seconds;
    };
    // CVodeReInit zeroes CVODE's counters, so they are added up before each restart.
    auto collect = [&] {
        long v = 0;
        CVodeGetNumSteps(mem, &v);               sol.steps += v;
        CVodeGetNumErrTestFails(mem, &v);        sol.failed_steps += v;
        CVodeGetNumNonlinSolvConvFails(mem, &v); sol.failed_steps += v;
        CVodeGetNumNonlinSolvIters(mem, &v);     sol.nonlin_iters += v;
    };
    /* Restarts keep the solver's t below one output interval: a stage time
       near t ~ 1e12 s is rounded to its ulp, which puts noise ~Delta eps t
       into exp(i Delta t), and the error estimate cannot get below it. The q
       equations carry no explicit t, so from t0 on A_x(t0 + s) = B_x(s)
       exp(i w_x t0), where B obeys the same dA_dt from s = 0 with B(0) = q(t0).
       The solver integrates B; the equations are untouched. */
    double t0 = 0.0, s = 0.0;                    // t = t0 + s, s the solver's t
    bool done = false;
    for (int step=0 ; step<opt.n_out && !done ; step++) {
        const double t_target = t_end * step / (opt.n_out - 1);
        if (s > 0.0) {
            collect();
            if (opt.slow) std::copy_n(to_q(s, State(as_complex(yd), as_complex(yd) + n)).begin(), n, as_complex(yd));
            CVodeReInit(mem, 0.0, y);            // GSL: gsl_odeiv2_driver_reset
            t0 = t;
            s = 0.0;
        }

        // advance from t to t_target, updating t and y in place; flag < 0 is a failure
        int flag = 0;
        double e = energy_now();
        while (t < t_target && flag == 0) {
            const double t_next = std::min(t_target, t + check);
            flag = CVode(mem, t_next - t0, y, &s, CV_NORMAL);
#if 0
            flag = gsl_odeiv2_driver_apply(drv, &s, t_next - t0, yd) == GSL_SUCCESS ? 0 : -1;
#endif
            t = t0 + s;
            e = energy_now();
            if (!std::isfinite(e)) { flag = -1; break; }
            if (t < t_target && crossed(e)) { done = true; break; }
            if (out_of_time()) { sol.out_of_time = done = true; break; }
        }
        if (flag < 0) { sol.failed = true; break; }
        sol.t.push_back(t);
        sol.y.push_back(opt.slow ? to_q(-t0, out) : out);   // A = B exp(i w t0)
        if (!done && crossed(e)) done = true;
    }

    sol.rhs_calls = params.rhs_calls;
    collect();
    CVodeFree(&mem);
    SUNNonlinSolFree(NLS);
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
