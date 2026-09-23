#pragma once
#include "angular.hpp"
#include "model.hpp"


constexpr std::array<int, 6> REFINE_LADDER = {1, 2, 4, 8, 16, 32};
constexpr double QUAD_TOL = 1e-4;

constexpr std::array<const char*, 8> GROUPS =
    {"A55", "A56", "A57", "A58", "A59", "A60", "A61", "A62"};

// Background on the merged integration grid.
struct Grid {
    eig::ArrayXd r, rho, P, Gamma_1, g, dg_dr, dlnrho_dlnr, dGamma1_dlnrho_s;
};

// One mode's radial functions on that same grid.
struct Fields {
    eig::ArrayXd xi_r, xi_h, div_xi, delta_Phi, ddelta_Phi_dr;
    double Lambda2, omega2;
};

using EF3   = std::array<const Eigenfunction*, 3>;
using F3    = std::array<Fields, 3>;
using Basis = std::array<double, 7>;        // T, F_a, F_b, F_c, G_a, G_b, G_c

struct KappaResult {
    double kappa, kappa_simpson;
    eig::ArrayXd r, dkappa_dr, cumulative;  // dkappa_dr already normalized by 2E_star
    std::array<double, 8> groups;           // Equations A55--A62
    int refine;
    bool converged;

    double quadrature_residual () const;
    // max |kappa(<r)| / |kappa(R)|. Above ~1e2 the integral is dominated by
    // cancellation and the result is not trustworthy.
    double cancellation () const;
    // Portion of |dkappa/dr| accumulated outside r > x_split R.
    double outer_fraction (double x_split = 0.8) const;
};

// refine = 0 runs REFINE_LADDER; any other value pins that single rung.
KappaResult kappa_abc (
    const Eigenfunction& a, const Eigenfunction& b, const Eigenfunction& c,
    std::array<int, 3> ms,
    int refine = 0, double tol = QUAD_TOL
);

// Every m of one triplet from one set of radial integrals. Returns kappa per
// entry of `ms`, in that order, and the refinement that satisfied `tol`.
std::pair<std::vector<double>, int>
kappa_all_m (
    const Eigenfunction& a, const Eigenfunction& b, const Eigenfunction& c,
    const std::vector<std::array<int, 3>>& ms,
    int refine = 0, double tol = QUAD_TOL
);

double kappa_from_basis (const Basis& b, const AngularFactors& ang);

std::pair<Grid, F3> on_grid (const EF3& efs, int refine);

// dkappa/dr split by which angular factor multiplies it. Every angular factor
// is a scalar, so kappa = sum_j <factor_j> I_j with the seven radial integrals
// independent of m.
std::array<eig::ArrayXd, 7> basis_densities (const Grid& star, const F3& f);

eig::ArrayXd A55 (const Grid&, const F3&, const AngularFactors&);
eig::ArrayXd A56 (const Grid&, const F3&, const AngularFactors&);
std::pair<eig::ArrayXd, eig::ArrayXd>
A56_parts ( // A56 split into its Lambda^2 xi_h and its -4 xi_r pieces.
    const Grid&, const F3&,
    const AngularFactors&
);
eig::ArrayXd A57 (const Grid&, const F3&, const AngularFactors&);
eig::ArrayXd A58 (const Grid&, const F3&, const AngularFactors&);
eig::ArrayXd A59 (const Grid&, const F3&, const AngularFactors&);
eig::ArrayXd A60 (const Grid&, const F3&, const AngularFactors&);
eig::ArrayXd A61 (const Grid&, const F3&, const AngularFactors&);
eig::ArrayXd A62 (const Grid&, const F3&, const AngularFactors&);
