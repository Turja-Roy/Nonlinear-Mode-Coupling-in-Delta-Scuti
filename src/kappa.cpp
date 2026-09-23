#include "kappa.hpp"
#include "numeric.hpp"

#include <algorithm>
#include <numeric>

namespace {

constexpr AngularFactors UNIT_T{1.0, 0, 0, 0, 0, 0, 0, 0};   // A55-A58, A62 carry only T

int cyc_jk (int i) { return (i + 1) % 3; }

eig::ArrayXd g4(const Grid& star) {
    return 4.0 * star.g + star.r * star.dg_dr; 
}

// Index one past x = 1. MESA grid continues into the atmosphere,
// but the integral stops at the photosphere
eig::Index cut_at_surface (const eig::ArrayXd& x) {
    return std::upper_bound(x.data(), x.data()+x.size(), 1.0) - x.data();
}

// refine-1 extra points inside every cell, endpoints preserved.
eig::ArrayXd subdivide (const eig::ArrayXd& x, int refine) {
    eig::ArrayXd out(refine * (x.size()-1) + 1);
    for (eig::Index i=0 ; i<x.size()-1 ; i++) {
        const double h = (x[i + 1] - x[i]) / refine;
        for (int f=0 ; f<refine ; f++)
            out[refine*i + f] = x[i] + h*f;
    }
    out[out.size()-1] = x[x.size()-1];
    return out;
}

using SM = eig::ArrayXd Star::*;
using GM = eig::ArrayXd Grid::*;
using EM = eig::ArrayXd Eigenfunction::*;
using FM = eig::ArrayXd Fields::*;

constexpr SM STAR_F[8] = {&Star::r, &Star::rho, &Star::P, &Star::Gamma_1, &Star::g,
                          &Star::dg_dr, &Star::dlnrho_dlnr, &Star::dGamma1_dlnrho_s};
constexpr GM GRID_F[8] = {&Grid::r, &Grid::rho, &Grid::P, &Grid::Gamma_1, &Grid::g,
                          &Grid::dg_dr, &Grid::dlnrho_dlnr, &Grid::dGamma1_dlnrho_s};
constexpr EM EIG_F[5]  = {&Eigenfunction::xi_r, &Eigenfunction::xi_h,
                          &Eigenfunction::div_xi, &Eigenfunction::delta_Phi,
                          &Eigenfunction::ddelta_Phi_dr};
constexpr FM FLD_F[5]  = {&Fields::xi_r, &Fields::xi_h, &Fields::div_xi,
                          &Fields::delta_Phi, &Fields::ddelta_Phi_dr};

}


/* ------------------
    Eqs. (A55)-(A62)
   ------------------
*/

eig::ArrayXd A55 (const Grid& star, const F3& f, const AngularFactors& ang) {
    const eig::ArrayXd bracket = star.Gamma_1 * (star.Gamma_1+1.0) + star.dGamma1_dlnrho_s;
    return ang.T * star.r.square() * star.P * bracket
           * f[0].div_xi * f[1].div_xi * f[2].div_xi;
}

eig::ArrayXd A56 (const Grid& star, const F3& f, const AngularFactors& ang) {
    eig::ArrayXd bracket = eig::ArrayXd::Zero(star.r.size());
    for (int i=0 ; i<3 ; i++) {
        const int j = cyc_jk(i), k = cyc_jk(i+1);
        bracket += f[j].div_xi * f[k].div_xi
               * (f[i].Lambda2 * f[i].xi_h - 4.0 * f[i].xi_r);
    }
    return ang.T * star.r * star.P * star.Gamma_1 * bracket;
}

std::pair<eig::ArrayXd, eig::ArrayXd>
A56_parts (const Grid& star, const F3& f, const AngularFactors& ang) {
    eig::ArrayXd hor = eig::ArrayXd::Zero(star.r.size());
    eig::ArrayXd rad = eig::ArrayXd::Zero(star.r.size());
    for (int i=0 ; i<3 ; i++) {
        const int j = cyc_jk(i), k = cyc_jk(i+1);
        const eig::ArrayXd pre = f[j].div_xi * f[k].div_xi;
        hor += pre * f[i].xi_h * f[i].Lambda2;
        rad += pre * -4.0 * f[i].xi_r;
    }
    const eig::ArrayXd scale = ang.T * star.r * star.P * star.Gamma_1;
    return {scale * hor, scale * rad};
}

eig::ArrayXd A57 (const Grid& star, const F3& f, const AngularFactors& ang) {
    const eig::ArrayXd drho_dlnr = star.rho * star.dlnrho_dlnr;
    return ang.T * drho_dlnr * g4(star) * f[0].xi_r * f[1].xi_r * f[2].xi_r;
}

eig::ArrayXd A58 (const Grid& star, const F3& f, const AngularFactors& ang) {
    eig::ArrayXd out = eig::ArrayXd::Zero(star.r.size());
    for (int i=0 ; i<3 ; i++)
        out += f[i].div_xi * f[cyc_jk(i)].xi_r * f[cyc_jk(i+1)].xi_r;
    return ang.T * star.rho * star.r * g4(star) * out;
}

eig::ArrayXd A59 (const Grid& star, const F3& f, const AngularFactors& ang) {
    const double w = f[0].omega2*ang.G_a + f[1].omega2*ang.G_b + f[2].omega2*ang.G_c;
    return -star.rho * star.r * f[0].xi_h * f[1].xi_h * f[2].xi_h * w;
}

eig::ArrayXd A60 (const Grid& star, const F3& f, const AngularFactors& ang) {
    const double F[3] = {ang.F_a, ang.F_b, ang.F_c};
    eig::ArrayXd out = eig::ArrayXd::Zero(star.r.size());
    for (int i=0 ; i<3 ; i++) {
        const int j = cyc_jk(i), k = cyc_jk(i+1);
        const double w = (f[i].omega2 - 3.0*f[j].omega2 - 3.0*f[k].omega2) * F[i]
                       - 2.0 * (f[j].omega2*F[j] + f[k].omega2*F[k]);
        out += f[i].xi_r * f[j].xi_h * f[k].xi_h * w;
    }
    return -star.rho * star.r * out;
}

eig::ArrayXd A61 (const Grid& star, const F3& f, const AngularFactors& ang) {
    const double F[3] = {ang.F_a, ang.F_b, ang.F_c};
    eig::ArrayXd out = eig::ArrayXd::Zero(star.r.size());
    for (int i=0 ; i<3 ; i++) {
        const int j = cyc_jk(i), k = cyc_jk(i+1);
        const double w = f[j].omega2*F[j] + f[k].omega2*F[k]
                       - 6.0*f[i].omega2*ang.T;
        out += f[i].xi_h * f[j].xi_r * f[k].xi_r * w;
    }
    return star.rho * star.r * out;
}

eig::ArrayXd A62 (const Grid& star, const F3& f, const AngularFactors& ang) {
    eig::ArrayXd out = eig::ArrayXd::Zero(star.r.size());
    for (int i=0 ; i<3 ; i++) {
        const int j = cyc_jk(i), k = cyc_jk(i+1);
        const eig::ArrayXd shell = star.dlnrho_dlnr * f[j].xi_r * f[k].xi_r
            + star.r * (f[j].xi_r*f[k].div_xi + f[k].xi_r*f[j].div_xi);
        out += shell * (star.r*f[i].ddelta_Phi_dr + 2.0*f[i].delta_Phi);
    }
    return ang.T * star.rho * out;
}


/* -----------------
    Basis densities
   -----------------
*/

std::array<eig::ArrayXd, 7> basis_densities (const Grid& star, const F3& f) {
    std::array<eig::ArrayXd, 3> P, S;
    for (int i=0 ; i<3 ; i++) {
        const int j = cyc_jk(i), k = cyc_jk(i+1);
        P[i] = star.rho * star.r * f[i].xi_r * f[j].xi_h * f[k].xi_h;   // A60's radial parts
        S[i] = star.rho * star.r * f[i].xi_h * f[j].xi_r * f[k].xi_r;   // A61's radial parts
    }
    const double w[3] = {f[0].omega2, f[1].omega2, f[2].omega2};

    std::array<eig::ArrayXd, 7> out;
    out[0] = A55(star, f, UNIT_T) + A56(star, f, UNIT_T) + A57(star, f, UNIT_T)
           + A58(star, f, UNIT_T) + A62(star, f, UNIT_T)
           - 6.0 * (w[0] * S[0] + w[1] * S[1] + w[2] * S[2]);       // A61's T remainder

    const eig::ArrayXd hhh = -star.rho * star.r * f[0].xi_h * f[1].xi_h * f[2].xi_h;
    for (int i=0 ; i<3 ; i++) {
        const int j = cyc_jk(i), k = cyc_jk(i+1);
        out[1 + i] = -((w[i] - 3.0 * w[j] - 3.0 * w[k]) * P[i]
                       - 2.0 * w[i] * (P[j] + P[k]))
                   + w[i] * (S[j] + S[k]);
        out[4 + i] = w[i] * hhh;
    }
    return out;
}

double kappa_from_basis (const Basis& b, const AngularFactors& ang) {
    const std::array<double, 7> a = ang.basis();
    double s = 0.0;
    for (int j=0 ; j<7 ; j++) s += b[j] * a[j];
    return s;
}


/* --------------------------
    One common spatial grid
   --------------------------
*/

std::pair<Grid, F3> on_grid (const EF3& efs, int refine) {
    Grid star;
    F3 f;
    for (int i=0 ; i<3 ; i++) {
        f[i].Lambda2 = efs[i]->Lambda2();
        f[i].omega2  = efs[i]->omega * efs[i]->omega;
    }

    const Star& s0 = *efs[0]->starptr;
    if (efs[0]->starptr == efs[1]->starptr && efs[1]->starptr == efs[2]->starptr) {
        // Shared grid: slice and use as they are. Exact, and the common case.
        const eig::Index n = cut_at_surface(s0.x);
        for (int q=0 ; q<8 ; q++) star.*GRID_F[q] = (s0.*STAR_F[q]).head(n);
        for (int i=0 ; i<3 ; i++)
            for (int q=0 ; q<5 ; q++)
                f[i].*FLD_F[q] = (*efs[i].*EIG_F[q]).head(n);

        if (refine > 1) {
            const eig::ArrayXd x = s0.x.head(n), xn = subdivide(x, refine);
            for (int q = 0; q < 8; ++q)
                star.*GRID_F[q] = num::Pchip(x, star.*GRID_F[q]).eval(xn);
            for (int i = 0; i < 3; ++i)
                for (int q = 0; q < 5; ++q)
                    f[i].*FLD_F[q] = num::CubicSpline(x, f[i].*FLD_F[q]).eval(xn);
        }
        return {star, f};
    }

    /* Distinct grids. The master grid is the union of the three, so every
       point has an exact (stellar) background value from at least one of them 
       and only the eigenfunctions are interpolated. A mode is smooth wherever GYRE
       declined to refine it, which is exactly where the other two add points. */
    std::array<eig::ArrayXd, 3> xs;
    eig::Index total = 0;
    for (int i = 0; i < 3; ++i) {
        const Star& s = *efs[i]->starptr;
        xs[i] = s.x.head(cut_at_surface(s.x));
        total += xs[i].size();
    }

    std::vector<double> u;
    u.reserve(total);
    for (int i=0 ; i<3 ; i++)
        u.insert(u.end(), xs[i].data(), xs[i].data() + xs[i].size());
    std::sort(u.begin(), u.end());
    u.erase(std::unique(u.begin(), u.end()), u.end());
    eig::ArrayXd x = eig::Map<const eig::ArrayXd>(u.data(), eig::Index(u.size()));
    if (refine > 1)
        x = subdivide(x, refine);

    // Concatenate the three structure samples, sort once, drop duplicate x.
    // Where a sample already holds a point the interpolant reproduces it, so
    // only genuinely new points are approximated.
    eig::ArrayXd x_all(total);
    for (int i=0,at=0 ; i<3 ; at+=xs[i].size(),i++)
        x_all.segment(at, xs[i].size()) = xs[i];

    std::vector<eig::Index> ord(total);
    std::iota(ord.begin(), ord.end(), eig::Index(0));
    std::stable_sort(ord.begin(), ord.end(),
                     [&](eig::Index p, eig::Index q) { return x_all[p] < x_all[q]; });
    std::vector<eig::Index> keep;
    keep.reserve(ord.size());
    for (size_t t = 0; t < ord.size(); ++t)
        if (t == 0 || x_all[ord[t]] > x_all[ord[t - 1]]) keep.push_back(ord[t]);

    eig::ArrayXd xk(eig::Index(keep.size())), vk(eig::Index(keep.size()));
    for (size_t t=0 ; t<keep.size() ; t++) xk[eig::Index(t)] = x_all[keep[t]];

    // PCHIP for the background because it has kinks; cubic for the
    // eigenfunctions because they do not. Not to swap them.
    eig::ArrayXd v_all(total);
    for (int q=0 ; q<8 ; q++) {
        for (int i=0,at=0 ; i<3 ; i++) {
            const Star& s = *efs[i]->starptr;
            v_all.segment(at, xs[i].size()) = (s.*STAR_F[q]).head(xs[i].size());
            at += xs[i].size();
        }
        for (size_t t=0 ; t<keep.size() ; t++) vk[eig::Index(t)] = v_all[keep[t]];
        star.*GRID_F[q] = num::Pchip(xk, vk).eval(x);
    }
    for (int i=0 ; i<3 ; i++)
        for (int q=0 ; q<5 ; q++) {
            const eig::ArrayXd v = (*efs[i].*EIG_F[q]).head(xs[i].size());
            f[i].*FLD_F[q] = num::CubicSpline(xs[i], v).eval(x);
        }
    return {star, f};
}


/* --------------------
    Quadrature ladder
   --------------------
*/

namespace {

// Spline and Simpson versions of the same seven integrals, already normalised.
std::pair<Basis, Basis> basis_pair (const Grid& star, const F3& f, double norm) {
    const std::array<eig::ArrayXd, 7> d = basis_densities(star, f);
    Basis spline, simp;
    for (int j=0 ; j<7 ; j++) {
        spline[j] = num::CubicSpline(star.r, d[j]).cumulative_integral()[star.r.size() - 1] * norm;
        simp[j]   = num::simpson(d[j], star.r) * norm;
    }
    return {spline, simp};
}

double rel (double a, double b) {
    return std::abs(a-b) / std::max(std::abs(a),1e-300); 
}

}  // namespace

std::pair<std::vector<double>, int>
kappa_all_m (const Eigenfunction& a, const Eigenfunction& b, const Eigenfunction& c,
            const std::vector<std::array<int, 3>>& ms, int refine, double tol) {
    const EF3 efs = {&a, &b, &c};
    const double norm = 1.0 / (2.0 * a.starptr->E_star());

    std::vector<AngularFactors> angs;
    angs.reserve(ms.size());
    for (const auto& m : ms) angs.push_back(angular_factors(a.l, b.l, c.l, m[0], m[1], m[2]));

    std::vector<double> out(ms.size());
    int n = refine;
    for (size_t rung = 0; rung < REFINE_LADDER.size(); ++rung) {
        n = refine ? refine : REFINE_LADDER[rung];
        const auto [star, f] = on_grid(efs, n);
        const auto [spline, simp] = basis_pair(star, f, norm);

        // Judged on kappa, not on the seven components: a near-null kappa can
        // come from cancellation between individually well-resolved pieces.
        double resid = 0.0;
        for (size_t q = 0; q < ms.size(); ++q) {
            out[q] = kappa_from_basis(spline, angs[q]);
            resid = std::max(resid, rel(out[q], kappa_from_basis(simp, angs[q])));
        }
        if (resid <= tol || refine) break;
    }
    return {out, n};
}

KappaResult kappa_abc (const Eigenfunction& a, const Eigenfunction& b,
                      const Eigenfunction& c, std::array<int, 3> ms,
                      int refine, double tol) {
    const EF3 efs = {&a, &b, &c};
    const double norm = 1.0 / (2.0 * a.starptr->E_star());
    const AngularFactors ang = angular_factors(a.l, b.l, c.l, ms[0], ms[1], ms[2]);

    Grid star;
    F3 f;
    Basis spline{};
    int n = refine;
    bool converged = refine != 0;
    for (size_t rung=0 ; rung<REFINE_LADDER.size() ; rung++) {
        n = refine ? refine : REFINE_LADDER[rung];
        std::tie(star, f) = on_grid(efs, n);
        Basis simp;
        std::tie(spline, simp) = basis_pair(star, f, norm);
        if (rel(kappa_from_basis(spline, ang), kappa_from_basis(simp, ang)) <= tol) {
            converged = true;
            break;
        }
        if (refine) break;
    }

    // dkappa/dr is the same density either way; the basis contraction is just
    // the group sum with the angular scalars pulled out.
    const std::array<eig::ArrayXd, 7> d = basis_densities(star, f);
    const std::array<double, 7> af = ang.basis();
    eig::ArrayXd total = eig::ArrayXd::Zero(star.r.size());
    for (int j = 0; j < 7; ++j) total += af[j] * d[j];
    total *= norm;

    const num::CubicSpline sp(star.r, total);
    const eig::ArrayXd cum = sp.cumulative_integral();

    KappaResult res;
    res.kappa         = cum[star.r.size() - 1];
    res.kappa_simpson = num::simpson(total, star.r);
    res.r             = star.r;
    res.dkappa_dr     = total;
    res.cumulative    = cum;
    res.refine        = n;
    res.converged     = converged;

    eig::ArrayXd (*const funcs[8])(const Grid&, const F3&, const AngularFactors&) =
        {A55, A56, A57, A58, A59, A60, A61, A62};
    for (int q=0 ; q<8 ; q++)
        res.groups[q] =
            num::CubicSpline(star.r, funcs[q](star, f, ang)).cumulative_integral()[star.r.size()-1] * norm;
    return res;
}

double KappaResult::quadrature_residual () const {
    return rel(kappa, kappa_simpson); 
}

double KappaResult::cancellation () const {
    return cumulative.abs().maxCoeff() / std::max(std::abs(kappa),1e-300);
}

double KappaResult::outer_fraction (double x_split) const {
    const eig::ArrayXd w = dkappa_dr.abs();
    const double total = num::trapz(w, r);
    if (total == 0.0) return 0.0;
    const double cutr = x_split * r[r.size() - 1];
    eig::Index lo = std::upper_bound(r.data(), r.data() + r.size(), cutr) - r.data();
    if (lo >= r.size() - 1) return 0.0;
    return num::trapz(w.tail(r.size() - lo), r.tail(r.size() - lo)) / total;
}
