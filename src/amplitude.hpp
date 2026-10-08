#pragma once
#include "kappa.hpp"
#include "triplets.hpp"

#include <complex>
#include <string>
#include <vector>

/* Amplitude equations for a network of resonant triplets. MW25's own equation
   (Reading/Notes/sections/theory-overview.tex, eq:amplitude) and its
   self-coupled special case, worked out explicitly in their four-mode network
   (Reading/Notes/sections/mw25-reading.tex, eq:mw25four):

       dq_a/dt = -(i w_a + g_a) q_a + 2 i w_a kappa conj(q_b) conj(q_c)   (b != c)

       dq_a/dt = -(i w_a + g_a) q_a +   i w_a kappa conj(q_b)^2          (b == c,
       dq_b/dt = -(i w_b + g_b) q_b + 2 i w_b kappa conj(q_a) conj(q_b)   "a -> b + b")

   Signed frequencies, no sum-slot: the triplet resonates when its three
   signed omegas sum to ~0, so exactly one of them (the sum mode, or its
   negation) carries the odd sign. A real mode x is two variables, x+ at +|w|
   and x- at -|w| with q_x- = conj(q_x+) (Schenk et al. 2002), and every
   resonance appears twice, (s-, p+, q+) and its mirror (s+, p-, q-); each
   variable is pushed by one of the two, so nothing is counted twice. Signing
   each mode once by its gamma, as MW25's examples do, works only for
   parametric and direct-sum triplets; a closure triplet d1 -> d2 + d3 then
   has |Delta| ~ 2w and does nothing. network_build writes the doubled set.
   The same three fields (a, b, c) work for every channel without relabelling.

   kappa is the bare coupling integral kappa_abc() returns. The explicit 2
   above is MW25's own combinatorial factor -- there are two orderings of two
   distinct partners -- and it collapses to 1 when the two partners are
   literally the same mode (b == c), because there is then only one ordering.
   Getting this wrong (always using 2) overcounts mode a's own equation by 2x
   in the self-coupled case; mw25-reading.tex's worked example is where these
   coefficients are checked against the paper.

   Slow amplitudes A
   -----------------
   q_x spins at the pulsation frequency, period of hours, while |q_x| changes
   over >= 100 d. Put q_x = A_x exp(-i w_x t) and integrate A instead: the
   spinning is carried by the exp, A holds the slow part, and a run that would
   need ~1e9 steps in q needs a few thousand in A. Substituting q_x = A_x
   exp(-i w_x t) into the q equation above (equivalently, MW25's own
   (epsilon, alpha) reduction, theory-overview.tex eq:slow, converted back to
   complex form) gives

       dA_a/dt = -g_a A_a + 2 i w_a kappa conj(A_b) conj(A_c) exp(i Delta t)   (b != c)

       dA_a/dt = -g_a A_a +   i w_a kappa conj(A_b)^2         exp(i Delta t)   (b == c)
       dA_b/dt = -g_b A_b + 2 i w_b kappa conj(A_a) conj(A_b) exp(i Delta t)

   with Delta = w_a + w_b + w_c (or w_a + 2 w_b when b == c) the same for
   every slot of the triplet -- unlike a sum-slot scheme, nothing here treats
   one mode's phase differently from another's. At exact resonance (Delta = 0)
   these carry no explicit t at all, so fixed points and thresholds can be
   read straight off them.

   Manley-Rowe is the check that catches a wrong sign or factor. With signed
   omega, N_x = |q_x|^2 / w_x carries the mode's own sign, and the equation's
   full a/b/c symmetry (unlike a sum-slot scheme, where a's equation looked
   different from b's and c's) makes the invariant a plain equality, not an
   alternating one. No damping, Delta = 0:

       b != c:  dN_a/dt = dN_b/dt = dN_c/dt
       b == c:  dN_b/dt = 2 dN_a/dt
*/

namespace amp {

using State = std::vector<std::complex<double>>;

struct Mode {
    std::string name;
    double omega;                                // signed, see above; not tied to gamma
    double gamma;                                // < 0 self-excited, > 0 damped
};

// Indices into Network::modes. b == c means a -> b + b (literally the same
// mode index, not two modes that merely share omega/gamma).
struct Triplet {
    int a, b, c;
    double kappa;
};

struct Options {
    int n_out = 2000;
    double rtol = 1e-12;
    double atol = 0.0;                           // 0 -> 1e-14 max|y0|, since |q| ~ 1e-6
    bool slow = true;                            // integrate A; false integrates q
    double e_max = 0.0;                          // > 0: stop when sum|y|^2 crosses it
    std::vector<int> frozen;                     // held at their initial amplitude
};

struct Solution {
    std::vector<double> t;
    std::vector<State> y;                        // A, or q when slow = false
};

class Network {
public:
    Network (std::vector<Mode> modes, std::vector<Triplet> triplets);

    int size () const { return int(modes_.size()); }
    const std::vector<Mode>& modes () const { return modes_; }
    const std::vector<Triplet>& triplets () const { return triplets_; }

    double detuning (const Triplet& tr) const;   // w_a + w_b + w_c

    // Right-hand sides of the equations above: given (t, q) they return dq/dt.
    State dq_dt (double t, const State& q) const;              // full q, fast phase kept
    State dA_dt (double t, const State& A) const;              // slow amplitudes A

    /* Same two, writing into a caller-supplied buffer. The solver rhs in
       integrate() owns the output buffer, which is the only reason these exist:
       the returning versions would allocate a State on every step, and a long
       run takes millions of them. */
    void dq_dt (double t, const State& q, State& out) const;
    void dA_dt (double t, const State& A, State& out) const;

    eig::ArrayXd energy (const State& y) const;  // E_x / E_star = |q_x|^2
    eig::ArrayXd action (const State& y) const;  // N_x = E_x / w_x
    State to_q (double t, const State& A) const;             // put the phase back

    /* e_max stops the run on the way up only, so a runaway shows as
       t.back() < t_end, the last output being the first past e_max; a network with no damped mode has no bounded state and
       without the cap the integrator just grinds against the blow-up.
       Freezing is only meaningful for the A equations. */
    Solution integrate (State y0, double t_end, Options opt = Options{}) const;

private:
    // phase = exp(i Delta t) for the A equations, 1 for the q equations, the
    // same value for every slot of the triplet.
    void add_triplet (const Triplet& tr, const State& y, State& out,
                      std::complex<double> phase) const;
    void add_triplets (double t, const State& y, State& out, bool slow) const;

    std::vector<Mode> modes_;
    std::vector<Triplet> triplets_;
    eig::ArrayXd omega_, gamma_;
    std::vector<double> delta_;                  // one per triplet, cached
};

// omega and gamma must already be signed (see Mode above).
Network three_mode (std::array<double, 3> omega, std::array<double, 3> gamma, double kappa,
                    std::array<std::string, 3> names = {"a", "b", "c"});

// Mode a decaying into two copies of one mode d: a -> d + d
Network self_coupled (double omega_a, double omega_d, double gamma_a, double gamma_d,
                      double kappa);

}  // namespace amp
