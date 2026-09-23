#pragma once
#include "angular.hpp"
#include "model.hpp"

#include <set>

/* Resonant triplet enumeration.

   omega is m-degenerate for a nonrotating star, so the detuning depends only
   on (l, n_pg). Enumeration works at that level and attaches the m
   combinations afterwards; kappa exploits the same split, computing each
   radial integral once per RadialTriplet and contracting it with the angular
   factors of every m combination.

   Sign classes: with all stored omega positive, (+,+,+) can never resonate
   here, and any two-minus assignment is the negation of a one-minus one, so
   |Delta| is unchanged. One mode carries s = -1 -- the sum mode, always the
   highest-frequency member.

   Slots a, b, c are bookkeeping, not roles: a is simply the sum mode. Which
   modes are parents and which are daughters is decided by linear stability
   alone -- see stab::channel() in stability.hpp. */

using Key = std::pair<int, int>;                 // (l, n_pg)

constexpr double DETUNING_CUT_DIMLESS = 0.15;    // in units of sqrt(GM/R^3)

struct RadialTriplet {
    Key sum_mode;
    std::array<Key, 2> pair;
    std::array<double, 3> omega;                 // signed, rad/s; sum mode first, negative
    double delta;                                // rad/s

    std::array<Key, 3> keys() const { return {sum_mode, pair[0], pair[1]}; }
    std::array<int, 3> ls() const;
    std::array<int, 3> ns() const;
    static constexpr std::array<int, 3> signs() { return {-1, +1, +1}; }

    std::vector<std::array<int, 3>> m_combinations() const;
    std::array<Mode, 3> modes(const ModeMap& efs, std::array<int, 3> ms) const;
};

/* Selection rules first, then |Delta| < cut, in rad/s.

   `sum_keys` restricts slot a, the sum mode, to those (l, n_pg); the pair
   still ranges over the whole net. That asymmetry is what makes a wide-net
   scan tractable: only a self-excited mode can pump a parametric pair, and
   the driven set is a few dozen modes against a few thousand. */
std::vector<RadialTriplet> enumerate_triplets(const ModeMap& efs, double cut,
                                              int l_max = 3,
                                              const std::set<Key>* sum_keys = nullptr);

long count_with_m(const std::vector<RadialTriplet>& triplets);

struct TripletSummary {
    long n_radial, n_with_m;
    double min_abs_delta_over_omega, median_abs_delta_over_omega;
    std::map<std::array<int, 3>, long> l_combinations;
};
TripletSummary summarise(const std::vector<RadialTriplet>& triplets);

// triplets_<tag>.csv, column for column.
struct TripletRow {
    int l_a, n_a, l_b, n_b, l_c, n_c;
    double omega_a, omega_b, omega_c, delta, abs_delta_over_omega_a;
    long n_m;                                    // -1 when m counting was skipped
};
std::vector<TripletRow> triplet_rows(const std::vector<RadialTriplet>& triplets,
                                     bool with_m = true);
