// GSL for ODE integration
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <limits>
#include <random>
#include <string>
#include <vector>

using cplx = std::complex<double>;
static const cplx I(0.0, 1.0);

struct Mode {
    std::string name;
    double omega;
    double gamma;
};

// Indices into the mode list. b == c means a -> b + b, one mode literally
// repeated (not two modes that merely share omega/gamma).
struct Triplet {
    int a, b, c;
    double kappa;
};

struct AmplEqParams {
    const std::vector<Mode>* modes;
    const std::vector<Triplet>* triplets;
};

double detuning (const std::vector<Mode>& modes, const Triplet& tr) {
    return modes[tr.a].omega + modes[tr.b].omega + modes[tr.c].omega;
}

/* -------------------------
    The Amplitude Equations
   ------------------------- */
int dq_dt (double t, const double y[], double dydt[], void* params) {
    (void)t;                                                // unused
    const auto* p = static_cast<AmplEqParams*>(params);     // pointer containing modes and triplets
    const std::vector<Mode>& modes = *p->modes;
    const std::vector<Triplet>& triplets = *p->triplets;
    const int n = int(modes.size());                        // number of modes

    const auto* q = reinterpret_cast<const cplx*>(y);       // input amplitudes
    auto* out = reinterpret_cast<cplx*>(dydt);              // output derivatives

    // Linear term: -(i omega + gamma) q
    // The linear term is separated because regardless of how many triplets we're integrating,
    // the linear term is always the same and only depends on the mode's own omega and gamma.
    //                               = - (i omega + gamma) q
    for (int m=0 ; m<n ; m++) out[m] = -(I * modes[m].omega + modes[m].gamma) * q[m];

    // Nonlinear term: 2 i omega kappa conj(q_j) conj(q_k), added once per
    // triplet the mode is in (e.g. 2, 10, 50 or more daughters -> that many
    // triplets, hence the loop).
    //
    // b == c (a -> b + b, one mode literally repeated) is not the general
    // case with the same index written twice: MW25's own worked example
    // (mw25-reading.tex eq:mw25four) shows the repeated mode's own equation
    // keeps the factor 2, but the *other* mode's equation drops to 1 -- there
    // is only one way to pick "both partners = b" instead of two orderings.
    // Silently using 2 for both would overcount mode a's driving by 2x.
    for (const Triplet& tr : triplets) {
        const int a=tr.a , b=tr.b , c=tr.c;
        if (b == c) {
            out[a] +=       I * modes[a].omega * tr.kappa * std::conj(q[b]) * std::conj(q[b]);
            out[b] += 2.0 * I * modes[b].omega * tr.kappa * std::conj(q[a]) * std::conj(q[b]);
            continue;
        }
        //     += 2 i omega kappa conj(q_j) conj(q_k)
        out[a] += 2.0 * I * modes[a].omega * tr.kappa * std::conj(q[b]) * std::conj(q[c]);
        out[b] += 2.0 * I * modes[b].omega * tr.kappa * std::conj(q[a]) * std::conj(q[c]);
        out[c] += 2.0 * I * modes[c].omega * tr.kappa * std::conj(q[a]) * std::conj(q[b]);
    }
    return GSL_SUCCESS;
}

// mu = |omega kappa| / sqrt(Delta^2 + gamma^2), q_driven = mu q_1 q_2.
double mu_of (double kappa, double omega, double delta, double gamma) {
    return std::abs(omega * kappa) / std::sqrt(delta*delta + gamma*gamma);
}

/* MW25 Eq. (6): energy at which a damped pair (b, c) pumped by a parent goes
   parametrically unstable. Needs both daughters actually damped. */
double threshold_energy (double kappa, double omega_b, double omega_c,
                         double gamma_b, double gamma_c, double delta) {
    if (gamma_b<=0.0 || gamma_c<=0.0 || kappa==0.0) return 0.0;

    const double g_sum = gamma_b + gamma_c;
    return gamma_b*gamma_c / (4.0 * omega_b*omega_c * kappa*kappa) *
           (1.0 + delta*delta / (g_sum*g_sum));
}

// MW25 Eq. (6)
double equilibrium_energy (double kappa, double omega_b, double omega_c,
                           double gamma_a, double gamma_b, double gamma_c, double delta) {
    if (gamma_b<=0.0 || gamma_c<=0.0 || kappa==0.0) return 0.0;

    const double g_sum = gamma_a + gamma_b + gamma_c;
    if (g_sum == 0.0) return std::numeric_limits<double>::infinity();
    return gamma_b*gamma_c / (4.0 * omega_b*omega_c * kappa*kappa) *
           (1.0 + delta*delta / (g_sum*g_sum));
}


static std::filesystem::path OUT;
static std::vector<std::array<std::string, 3>> LINES;        // panel, label, value

// Write a line to lines.csv, which is sorted and written at the end of main().
void line (const std::string& panel, const std::string& label, double v) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    LINES.push_back({panel, label, buf});
}

/* --------------------
    The ODE Integrator
   -------------------- */
/* Integrate one network with GSL's adaptive Prince-Dormand 8(9) stepper and
   write t plus one energy column per mode to <out>/<name>.csv. e_max > 0
   stops the run the first time sum|q|^2 crosses it (a runaway then shows as
   fewer rows than n_out, instead of the integrator grinding against the
   blow-up); t_scale divides the written time column (days for Figure 6). */
void run_panel (const std::string& name, const std::vector<Mode>& modes,
                const std::vector<Triplet>& triplets, std::vector<cplx> q0,
                double t_end, int n_out, double rtol, double e_max = 0.0,
                double t_scale = 1.0) {
    const int n = int(modes.size());        // Number of modes, also the size of q0.
    double qmax = std::max_element(q0.begin(), q0.end(),
        [] (const cplx& a, const cplx& b) { return std::abs(a) < std::abs(b); })->real();
    // Find largest initial amplitude to set the absolute tolerance relative to it.
    const double atol = 1e-14 * qmax;

    // Just packages the pointers to pass to GSL
    AmplEqParams params{&modes, &triplets};

    // System description for GSL: a struct holding all the ingredients together
    // dq_dt is the function, nullptr means no Jacobian, 2*n because real + imaginary parts, and params is the pointer to the parameters.
    gsl_odeiv2_system sys{dq_dt, nullptr, size_t(2 * n), &params};

    // The GSL ODE driver: onws the stepper and current state
    // Algorithm: Prince-Dormand 8(9), Runge-Kutta with adaptive step size, 8th order with 9th order error estimate.
    // Initial step: t_end / (10*n_out), Tolerances: atol (absolute) and rtol (relative)
    // gsl_odeiv2_driver* drv = gsl_odeiv2_driver_alloc_y_new(
    //     &sys, gsl_odeiv2_step_rk8pd, t_end / (10.0 * n_out), atol, rtol);
    // gsl_odeiv2_driver* drv = gsl_odeiv2_driver_alloc_y_new(
    //     &sys, gsl_odeiv2_step_rk4imp, t_end / (10.0 * n_out), atol, rtol);
    gsl_odeiv2_driver* drv = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_msadams, t_end / (10.0 * n_out), atol, rtol);

    // GSL works with raw doubles, y[2n] = {Re(q0[0]), Im(q0[0]), Re(q0[1]), Im(q0[1]), ...}
    std::vector<double> y(size_t(2*n));
    std::memcpy(y.data(), q0.data(), size_t(2*n) * sizeof(double)); // copy initial conditions

    // File I/O:
    std::FILE* f = std::fopen((OUT / (name + ".csv")).c_str(), "w");
    std::fprintf(f, "t");
    for (const Mode& m : modes) std::fprintf(f, ",E_%s", m.name.c_str());
    std::fprintf(f, "\n");

    // The integration loop: t goes from 0 to t_end, n_out times.
    double t = 0.0, e_prev = -1.0, t_last = 0.0;
    size_t rows = 0;
    bool stopped = false;
    for (int step=0 ; step<n_out ; step++) {
        const double t_target = t_end*step / (n_out-1);

        // GSL's driver_apply integrates from t to t_target, updating y in place.
        bool success = gsl_odeiv2_driver_apply(drv, &t, t_target, y.data());

        if (t_target>t && success != GSL_SUCCESS) {
            stopped = true; break;
        }
        if (e_max > 0.0) {
            // Total energy: |q|^2 = sum |q_i|^2 = sum (Re(q_i)^2 + Im(q_i)^2)
            double e = std::accumulate(y.begin(), y.end(), 0.0,
                                       [](double acc, double v) { return acc + v*v; });
            if (e_prev>=0.0 && e_prev<=e_max && e>e_max) {
                stopped = true; break; 
            }
            e_prev = e;
        }
        std::fprintf(f, "%.17g", t / t_scale);
        for (int m=0 ; m<n ; m++)
            std::fprintf(f, ",%.17g", y[2 * m] * y[2 * m] + y[2 * m + 1] * y[2 * m + 1]);
        std::fprintf(f, "\n");
        ++rows;
        t_last = t;
    }
    std::fclose(f);
    gsl_odeiv2_driver_free(drv);
    std::printf("  %-10s %5zu rows, t to %.4g%s\n", name.c_str(), rows, t_last / t_scale,
               stopped ? "  (stopped early: e_max)" : "");
}


/* ---------------
       Figures
   ---------------*/
void fig1 () {
    // (a) two identical parents a, b driving daughter c;
    const double KAPPA = 1.000;
    {
        // Values that matches closest to the paper's figures:
        // const double omega_a = -1.000, gamma_a = -0.010;
        // const double omega_b = -1.000, gamma_b = -0.010;
        // const double omega_c = 2.000, gamma_c = 0.100;

        // Values mentioned in the figure caption:
        const double omega_a = -0.5, gamma_a = -0.001;
        const double omega_b = -0.5, gamma_b = -0.001;
        const double omega_c = 1.0, gamma_c = 0.01;

        const std::vector<Mode> modes = {
            {"c", omega_c, gamma_c}, {"a", omega_a, gamma_a}, {"b", omega_b, gamma_b}};
        const std::vector<Triplet> triplets = {{0, 1, 2, KAPPA}};
        run_panel("fig1a", modes, triplets, {1e-3, 1e-3, 1e-3}, 6000.0, 8000, 1e-11);
    }
    // (b) two distinct parents a, b: the extra degree of freedom breaks it.
    {
        // Values that matches closest to the paper's figures:
        // const double omega_a = -1.000, gamma_a = -0.010;
        // const double omega_b = -1.004, gamma_b = -0.010;
        // const double omega_c = 2.000, gamma_c = 0.100;

        // Values mentioned in the figure caption:
        const double omega_a = -0.500, gamma_a = -0.0011;
        const double omega_b = -0.502, gamma_b = -0.001;
        const double omega_c = 1.000, gamma_c = 0.01;

        const std::vector<Mode> modes = {
            {"c", omega_c, gamma_c}, {"a", omega_a, gamma_a}, {"b", omega_b, gamma_b}};
        const std::vector<Triplet> triplets = {{0, 1, 2, KAPPA}};
        run_panel("fig1b", modes, triplets, {1e-3, 1e-3, 1e-3}, 6000.0, 8000, 1e-11, 1e2);
    }
}

// Two parents, N damped daughters, each an independent direct triplet.
void fig2 () {
    std::mt19937 gen(std::random_device{}());
    auto within_1pct = [&] (double b) { 
        return b * std::uniform_real_distribution<double>(0.99, 1.01)(gen); 
    };

    const double omega_a = -0.5, gamma_a = -0.001;
    const double omega_b = -0.5, gamma_b = -0.0011;
    const double W_D = 1.0, G_D = 0.01, K_D = 1.0;
    for (int N : {2, 10, 50}) {
        std::vector<Mode> modes = {{"a", omega_a, gamma_a}, {"b", omega_b, gamma_b}};
        std::vector<Triplet> triplets;
        for (int i=0 ; i<N ; i++) {
            const double wd = within_1pct(W_D), gd = within_1pct(G_D), kd = within_1pct(K_D);

            modes.push_back({"d" + std::to_string(i), wd, gd});
            triplets.push_back({2 + i, 0, 1, kd});
        }
        run_panel("fig2_N" + std::to_string(N), modes, triplets,
                 std::vector<cplx>(size_t(N + 2), 1e-7), 15000.0, 9000, 1e-10, 1e0);
    }
}

/* Parents a, b. Direct daughter c (a + b -> c). Self-coupled parametric
   daughter d (a -> d + d and b -> d + d), one mode literally repeated in each
   parametric triplet -- MW25's own four-mode network (mw25-reading.tex
   eq:mw25four), not two modes that merely share omega/gamma. */
const double omega[4] = {-1.000, -1.001, 2.0, 0.5};         // a, b, c, d
const double GAMMA_DEFAULT[4] = {-0.01, -0.01, 0.1, 0.1};   // a, b, c, d

void mixed (double k_direct, double k_param, const double g[4], std::vector<Mode>& modes,
            std::vector<Triplet>& triplets) {
    modes = {{"a", omega[0], g[0]}, {"b", omega[1], g[1]}, {"c", omega[2], g[2]},
             {"d", omega[3], g[3]}};
    triplets = {{2, 0, 1, k_direct},         // a + b -> c
                {0, 3, 3, k_param},          // a -> d + d
                {1, 3, 3, k_param}};         // b -> d + d
}

// E_a_th/E_a_eq: MW25 Eq. (6)/A7 for the parent a, pumped by its self-coupled
// daughter pair (d, d) -- not the daughter's own threshold, the parent's.
void mixed_lines (const std::string& tag, const std::vector<Mode>& modes,
                  const std::vector<Triplet>& triplets, double k_param, const double g[4]) {
    const double Delta = detuning(modes, triplets[1]);
    line(tag, "E_a_th", threshold_energy(k_param, omega[3], omega[3], g[3], g[3], Delta));
    line(tag, "E_a_eq", equilibrium_energy(k_param, omega[3], omega[3], g[0], g[3], g[3], Delta));
}

void fig3 () {
    std::vector<Mode> modes;
    std::vector<Triplet> triplets;
    mixed(1.0, 1.0, GAMMA_DEFAULT, modes, triplets);
    mixed_lines("3", modes, triplets, 1.0, GAMMA_DEFAULT);
    run_panel("fig3", modes, triplets, {1e-3, 1e-3, 1e-3, 1e-3}, 2000.0, 20000, 1e-11);
}

void fig4 () {
    // (a) parametric coupling 100x stronger: E_th scales as kappa^-2
    {
        std::vector<Mode> modes;
        std::vector<Triplet> triplets;
        mixed(1.0, 100.0, GAMMA_DEFAULT, modes, triplets);
        mixed_lines("4a", modes, triplets, 100.0, GAMMA_DEFAULT);
        run_panel("fig4a", modes, triplets, {3e-7, 3e-7, 3e-10, 3e-7}, 5000.0, 20000, 1e-11);
    }
    // (b) every linear rate 100x smaller: the cycle stretches out
    {
        double g[4];
        for (int i=0 ; i<4 ; i++) g[i] = 0.01 * GAMMA_DEFAULT[i];
        std::vector<Mode> modes;
        std::vector<Triplet> triplets;
        mixed(0.5, 1.0, g, modes, triplets);
        mixed_lines("4b", modes, triplets, 1.0, g);
        run_panel("fig4b", modes, triplets, {3e-7, 3e-7, 3e-10, 3e-7}, 800000.0, 800000, 1e-11);
    }
}

void fig5 () {
    // (a) E_c = |q_c|^2 and \mu^2 E_a E_b using fig3 parameters
    {
        std::vector<Mode> modes;
        std::vector<Triplet> triplets;
        mixed(1.0, 1.0, GAMMA_DEFAULT, modes, triplets);
        line("5a", "mu", mu_of(1.0, omega[2], detuning(modes, triplets[0]), GAMMA_DEFAULT[2]));
        run_panel("fig5a", modes, triplets, {1e-3, 1e-3, 1e-3, 1e-3}, 2000.0, 20000, 1e-11);
    }
    // (a) Same functions but with 30% of the driving/damping rates for the parents
    {
        double g[4] = {0.3*GAMMA_DEFAULT[0], 0.3*GAMMA_DEFAULT[1], GAMMA_DEFAULT[2], GAMMA_DEFAULT[3]};
        std::vector<Mode> modes;
        std::vector<Triplet> triplets;
        mixed(1.0, 1.0, g, modes, triplets);
        line("5b", "mu", mu_of(1.0, omega[2], detuning(modes, triplets[0]), GAMMA_DEFAULT[2]));
        run_panel("fig5b", modes, triplets, {1e-2, 1e-2, 1e-2, 1e-2}, 40000.0, 95000, 1e-11);
    }
}

void fig6() {
    constexpr double W_A = 1.17e-3, G_A = -4.8e-9, G_B = -1.1e-8;
    constexpr double W_D1 = 1.09e-3;
    constexpr double KP[2] = {0.85, 3.85}, DP[2] = {9.8e-8, 1.0e-4};
    (void)W_D1; (void)KP;     // fig6b scaffolding; re-enable when 5-mode panel returns
    const double KD = 6.3, DD = 7.6e-5, G_C = 2.8e-5;
    // const double KD  = v == 3 ? 63.0 : 6.3;
    // const double DD  = v == 1 || v == 2 ? 7.6e-6 : 7.6e-5;
    // const double G_C = v == 2 ? 2.8e-6 : 2.8e-5;

    // Magnitudes only; W_A here plays the same role as |omega_a| always did,
    // Delta = w_b+w_d1+w_d2 (say) reducing to the same arithmetic whether
    // omega_a is carried as a signed value or as a magnitude subtracted off.
    const double w_d2 = W_A + DP[0] - W_D1;
    const double w_b  = W_D1 + w_d2 + DP[1];
    const double w_c  = W_A + w_b - DD;

    // a, b are the self-excited parents (gamma < 0): signed negative here,
    // same magnitudes as above. c, d1, d2 are damped daughters: positive.
    // Seeded below the direct-coupling saddle (E ~ 5e-9) so the parent with
    // the larger omega and driving, b, is the one that runs away (MW25 fn. 1).
    run_panel("fig6a", {{"a", -W_A, G_A}, {"b", -w_b, G_B}, {"c", w_c, G_C}},
              {Triplet{2, 0, 1, KD}}, {1e-4, 1e-4, 3e-6}, 1e10, 3000, 1e-10, 1e0);
    // panel("fig6a", Network({{"a", -W_A, G_A}, {"b", -w_b, G_B}, {"c", w_c, G_C}},
    //                        {Triplet{2, 0, 1, KD}}),
          // State{1e-4, 1e-4, 3e-6}, t_end, 3000, 1e-10, 1e0, {}, DAY);

    // const Network net({{"a", -W_A, G_A}, {"b", -w_b, G_B}, {"c", w_c, G_C},
    //                    {"d1", W_D1, G_D1}, {"d2", w_d2, G_D2}},
    //                   {Triplet{2, 0, 1, KD},
    //                    Triplet{0, 3, 4, KP[0]},
    //                    Triplet{1, 3, 4, KP[1]}});
    // for (int i = 0; i < 2; ++i)
    //     line("6b", i ? "E_th_b" : "E_th_a",
    //          stab::threshold_energy(KP[i], W_D1, w_d2, G_D1, G_D2, DP[i]));
    // panel("fig6b", net, State{1e-4, 1e-4, 3e-6, 1e-5, 1e-5}, t_end, 4000, 1e-8, 1e0, {}, DAY);
    // std::printf("  variant %d: mu = %.0f\n", v, stab::mu(KD, w_c, DD, G_C));
}


int main (int argc, char** argv) {
    OUT = argc > 1 ? argv[1] : "out/mw25_original";
    std::filesystem::create_directories(OUT);
    std::printf("MW25 figures (original equations) -> %s\n", OUT.c_str());

    fig1();
    fig2();
    fig3();
    fig4();
    fig5();
    fig6();

    std::sort(LINES.begin(), LINES.end());
    std::FILE* f = std::fopen((OUT / "lines.csv").c_str(), "w");
    std::fprintf(f, "panel,label,value\n");
    for (const auto& l : LINES) std::fprintf(f, "%s,%s,%s\n", l[0].c_str(), l[1].c_str(), l[2].c_str());
    std::fclose(f);
    std::printf("  lines.csv  %zu entries\n", LINES.size());
    return 0;
}
