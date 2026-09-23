#pragma once
#include "kappa.hpp"
#include "triplets.hpp"

#include <complex>
#include <string>
#include <vector>

/* Amplitude equations for a network of resonant triplets

    a -> b + c
    a = sum slot
    b + c -> a
    a = sum slot

   Notes on conventions:
   ---------------------
   Paper uses a convention where parents have negative frequencies and children positive.
   In our code, a triplet is stored as (a, b, c) with all freqs being positive and a "sum mode".
   a being the sum mode (omega_a ~= omega_b + omega_c). So, in parametric coupling (a -> b + c),
   a is the parent while b and c are children, and in direct coupling (b + c -> a), a is the child.
   Detuning of that triplet:

       Delta = omega_b + omega_c - omega_a.

   The modes evolve as q_x(t) = A_x(t) exp(-i omega_x t). So, if we have an omega_a < 0,
   the conjugate is q_a(t) = A_a(t) exp(+i |omega_a| t).
   ** The conjugate of a mode is the same as the mode with negative frequency. **

   Hence, our amplitude equations become (with all positive frequencies):
   ----------------------------------------------------------------------
       dq_a/dt = -(i w_a + g_a) q_a + 2 i w_a kappa q_b    q_c
       dq_b/dt = -(i w_b + g_b) q_b + 2 i w_b kappa q_a  conj(q_c)
       dq_c/dt = -(i w_c + g_c) q_c + 2 i w_c kappa q_a  conj(q_b)

   Where the conjugates come from: a drive must oscillate like the mode it
   drives. With the exp(-i w t) convention,

       q_b q_c        ~ exp[-i(w_b + w_c)t] ~ exp(-i w_a t)   -> drives a
       q_a conj(q_c)  ~ exp[-i(w_a - w_c)t] ~ exp(-i w_b t)   -> drives b

   so the sum mode is driven by a plain product, and a pair mode is driven by
   the sum mode times the conjugate of its partner. Conjugating the partner is
   exactly what turns a sum frequency into a difference frequency. kappa is
   real here, and kappa* is written only to keep the convention visible.

   Slow amplitudes A
   -----------------
   q_x spins at the pulsation frequency, period of hours, while |q_x| changes
   over >= 100 d. Put q_x = A_x exp(-i w_x t) and integrate A instead: the
   spinning is carried by the exp, A holds the slow part, and a run that would
   need ~1e9 steps in q needs a few thousand in A. Same physics, phase factored
   out. The -i w_x q_x term cancels the fast phase of dq_x/dt, and each product
   is left with the slow phase exp(-+ i Delta t):

       dA_a/dt = -g_a A_a + 2 i w_a kappa A_b   A_c        exp(-i Delta t)
       dA_b/dt = -g_b A_b + 2 i w_b kappa A_a conj(A_c)    exp(+i Delta t)
       dA_c/dt = -g_c A_c + 2 i w_c kappa A_a conj(A_b)    exp(+i Delta t)

   At exact resonance (Delta = 0) these carry no explicit t at all, so fixed
   points and thresholds can be read straight off them.

   The factor 2 and the repeated-mode case
   ---------------------------------------
   The 2 counts the two orderings of the two partners. When the partners are
   the same mode it becomes 1, and when b and c are the same mode d
   (a -> d + d) the b and c equations are ONE equation, not two:

       dA_a/dt = -g_a A_a +   i w_a kappa  A_d^2        exp(-i Delta t)
       dA_d/dt = -g_d A_d + 2 i w_d kappa A_a conj(A_d) exp(+i Delta t)

   Adding a b-equation and a c-equation separately here double counts and is
   the one mistake this file exists to make impossible to overlook.

   Manley-Rowe fixes those factors and is the check that catches a wrong sign
   or factor. With N_x = |q_x|^2 / w_x, no damping and Delta = 0:

       distinct:  dN_a/dt = -dN_b/dt = -dN_c/dt
       repeated:  2 dN_a/dt = -dN_d/dt
*/

namespace amp {

using State = std::vector<std::complex<double>>;

struct Mode {
    std::string name;
    double omega;                                // > 0
    double gamma;                                // < 0 self-excited, > 0 damped
};

// Indices into Network::modes. b == c means a -> d + d.
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

    double detuning (const Triplet& tr) const;   // w_b + w_c - w_a

    // Right-hand sides of the equations above: given (t, q) they return dq/dt.
    State dq_dt (double t, const State& q) const;              // full q, fast phase kept
    State dA_dt (double t, const State& A) const;              // slow amplitudes A

    /* Same two, writing into a caller-supplied buffer. odeint calls its system
       as sys(y, dydt, t) and owns dydt, which is the only reason these exist:
       the returning versions would allocate a State on every step, and a long
       run takes millions of them. */
    void dq_dt (double t, const State& q, State& out) const;
    void dA_dt (double t, const State& A, State& out) const;

    eig::ArrayXd energy (const State& y) const;  // E_x / E_star = |q_x|^2
    eig::ArrayXd action (const State& y) const;  // N_x = E_x / w_x
    State to_q (double t, const State& A) const;             // put the phase back

    /* e_max stops the run on the way up only, so a runaway shows as
       t.back() < t_end; a network with no damped mode has no bounded state and
       without the cap the integrator just grinds against the blow-up.
       Freezing is only meaningful for the A equations. */
    Solution integrate (State y0, double t_end, Options opt = Options{}) const;

private:
    // phase = exp(-i Delta t) for the A equations, 1 for the q equations.
    void add_triplet (const Triplet& tr, const State& y, State& out,
                      std::complex<double> phase) const;

    std::vector<Mode> modes_;
    std::vector<Triplet> triplets_;
    eig::ArrayXd omega_, gamma_;
    std::vector<double> delta_;                  // one per triplet, cached
};

Network three_mode (std::array<double, 3> omega, std::array<double, 3> gamma, double kappa,
                    std::array<std::string, 3> names = {"a", "b", "c"});

// Sum mode a decaying into two copies of one mode d: a -> d + d
Network self_coupled (double omega_a, double omega_d, double gamma_a, double gamma_d,
                      double kappa);

/* Network over the union of several RadialTriplets. A mode shared between
   triplets is one node whatever role it plays in each.
   ms_list empty picks, per triplet, the m combination of largest |kappa|. */
Network from_triplets (const std::vector<RadialTriplet>& triplets, const ModeMap& efs,
                       const std::vector<std::array<int, 3>>& ms_list = {},
                       const std::map<Key, double>& gamma_override = {});

}  // namespace amp
