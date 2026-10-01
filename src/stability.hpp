#pragma once
#include "kappa.hpp"
#include "triplets.hpp"

#include <filesystem>

/* What one triplet does at linear order: who drives whom, how fast, and at
   what energy the nonlinear coupling turns on.

   gamma < 0 is a parent (self-excited), gamma > 0 a daughter (damped). Slot a
   is the sum mode, w_a ~= w_b + w_c; it is a parent in some channels and the
   daughter in others, so the roles come from the gammas, never from the slot.

   Energies below are E / E_star unless an E_star is passed in, in which case
   they come back in erg. */

namespace stab {

constexpr double NONADIABATIC_CYCLES = 2.0 * M_PI;

// Steady response of a driven mode x to a beat at detuning Delta:
//     mu = |w_x kappa| / sqrt(Delta^2 + gamma_x^2),   q_x = mu q_y q_z.
double mu (double kappa, double omega, double delta, double gamma);

// Growth rate of a daughter pair pumped by a parent of amplitude |q_a|:
//     Gamma = 2 |kappa| sqrt(|w_b w_c|) |q_a|.
double parametric_growth_rate (double kappa, double omega_b, double omega_c, double q_a);

/* Parent energy at which the daughter pair (b, c) goes parametrically unstable:

       E_th = gamma_b gamma_c / (4 kappa^2 |w_b w_c|) [1 + Delta^2/(gamma_b+gamma_c)^2].

   Zero when either daughter is itself driven (gamma <= 0): the pair is then
   linearly unstable and no threshold applies. */
double threshold_energy (double kappa, double omega_b, double omega_c,
                         double gamma_b, double gamma_c, double delta,
                         double E_star = 1.0);

/* Where the parent settles, MW25 eq. A7 -- E_th with the pair damping in the
   detuning penalty replaced by the triplet total:

       E_eq = gamma_b gamma_c / (4 kappa^2 |w_b w_c|)
              [1 + Delta^2/(gamma_a+gamma_b+gamma_c)^2],

   which is why a driven mode saturates near the amplitude at which it first
   went nonlinearly unstable. Infinite if the three gammas cancel: with no net
   dissipation the phase cannot lock and there is nothing to settle onto. */
double equilibrium_energy (double kappa, double omega_b, double omega_c,
                           double gamma_a, double gamma_b, double gamma_c,
                           double delta, double E_star = 1.0);

/* Largest E_th any (gamma_b, gamma_c) can give at |Delta| >> gamma: the
   bracket is maximised at gamma_b = gamma_c, leaving

       E_th,max = Delta^2 / (16 kappa^2 |w_b w_c|).

   An upper bound, which is what to quote while gamma is the least certain
   input. */
double threshold_energy_ceiling (double kappa, double omega_b, double omega_c,
                                 double delta, double E_star = 1.0);

/* All three slot energies at the A7 fixed point. Writing A7 cyclically and
   dividing any two collapses to

       E_a gamma_a / w_a = E_b gamma_b / w_b = E_c gamma_c / w_c,

   so E_a plus the ratios w_i/gamma_i gives the set -- and doing it that way
   avoids the sign trap that a term-by-term expansion falls into (for a
   parametric triplet gamma_a gamma_c and w_a w_c are both negative). NaN when
   the three w_i/gamma_i do not share a sign: the fixed point then puts a
   negative energy somewhere and there is nothing to settle onto. */
std::array<double, 3> equilibrium_energies (double kappa, std::array<double, 3> omega,
                                            std::array<double, 3> gamma,
                                            double delta, double E_star = 1.0);

// Shell where w t_thermal < 2 pi, i.e. where an adiabatic eigenfunction is no
// longer justified -- r/R >~ 0.96 for a delta Sct.
eig::Array<bool, eig::Dynamic, 1> nonadiabatic_shell (const Star& s, double omega);

// Share of |dkappa/dr| accumulated inside that shell: a physics caveat on the
// triplet, not a numerical one.
double nonadiabatic_fraction (const KappaResult& res, const Star& s, double omega);

/* Which channel the triplet can run, from the signs of the gammas alone:

     parametric   a is the only parent, b and c are its daughters
     direct-sum   b and c are parents, driving daughter a at w_b + w_c
     direct-diff  a and one of b, c are parents, driving the other at |w|
     all-driven   three parents, no daughter to absorb the energy
     all-damped   no parent
     inactive     a lone parent among b, c; its daughters are not in this triplet */
enum class Channel {
    parametric, direct_sum, direct_diff, all_driven, inactive, all_damped
};
constexpr std::array<Channel, 6> CHANNELS = {
    Channel::parametric, Channel::direct_sum, Channel::direct_diff,
    Channel::all_driven, Channel::inactive,   Channel::all_damped};

const char* channel_name (Channel c);
Channel channel (double gamma_a, double gamma_b, double gamma_c);

// Slot (0, 1, 2) of the daughter a direct channel drives; -1 if none runs.
// Not the same question as argmax(mu) -- the two disagree on ~13% of rows.
int daughter_slot (double gamma_a, double gamma_b, double gamma_c);

struct TripletObservables {
    std::array<int, 3> ms;
    double kappa;
    std::array<double, 3> mu;                    // per driven mode
    double E_threshold, E_equilibrium;           // erg
    int refine;

    double mu_max () const { return std::max({mu[0], mu[1], mu[2]}); }
};

// m000 keeps only m = (0, 0, 0). The radial integrals are m-independent, so
// the whole m set costs one 3j evaluation per combination.
std::vector<TripletObservables> observables (const RadialTriplet& t, const ModeMap& efs,
                                             bool m000 = false);

// One row of observables_<tag>.csv.
struct Row {
    int l_a, n_a, m_a, l_b, n_b, m_b, l_c, n_c, m_c;
    double omega_a, omega_b, omega_c;
    double gamma_a, gamma_b, gamma_c;
    double delta, kappa;
    double mu_a, mu_b, mu_c, mu_max;
    double E_th_over_E_star, E_eq_over_E_star, E_th_ceiling_over_E_star;
    bool detuning_dominated;
    int refine;
};

// One row per (triplet, m combination), ranked by mu_max descending.
// jobs > 1 spreads the kappa quadrature over OpenMP threads.
std::vector<Row> build_rows (const std::vector<RadialTriplet>& triplets, const ModeMap& efs,
                             const Star& s, bool m000 = false, int jobs = 1);

/* kappa at m = (0, 0, 0) per unordered (l, n) triple. The sign assignment
   does not enter kappa (only omega^2 does), so one value serves every channel
   that uses the triple. Legs missing from `kap` are integrated on `jobs`
   threads and appended to the TSV `cache` as they finish, so a crash loses
   only the integrals in flight. Delete the cache if the GYRE details change. */
using Leg = std::array<Key, 3>;                  // sorted
using KappaCache = std::map<Leg, std::pair<double, int>>;    // kappa, refine
Leg leg_of (Key a, Key b, Key c);
KappaCache load_kappa_cache (const std::filesystem::path& p);
void kappa_m000 (const std::vector<Leg>& legs, const ModeMap& efs,
                 const std::filesystem::path& cache, KappaCache& kap, int jobs);

// Column names and the values that go under them, kept side by side so they
// cannot drift apart. Use these rather than writing the row by hand.
extern const std::array<const char*, 26> ROW_COLUMNS;
void write_rows (const std::filesystem::path& p, const std::vector<Row>& rows,
                 bool force = true);

}  // namespace stab
