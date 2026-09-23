// g++ -std=c++17 -I/usr/include/eigen3 -Isrc tests/test_kappa.cpp
//     src/kappa.cpp src/angular.cpp src/numeric.cpp -lgsl -lgslcblas -o t && ./t
#include "kappa.hpp"

#include <cassert>
#include <set>
#include <cmath>
#include <cstdio>

static void close(double a, double b, double tol, const char* what) {
    if (std::abs(a - b) > tol * (1.0 + std::abs(b))) {
        std::fprintf(stderr, "FAIL %s: %.17g != %.17g\n", what, a, b);
        assert(false);
    }
}

// Smooth analytic stand-in for a star. The identities under test are algebraic,
// so the fields only have to be nonzero and differentiable.
static std::shared_ptr<Star> fake_star(eig::Index n) {
    auto s = std::make_shared<Star>();
    s->M = 4.0e33;
    s->R = 1.4e11;
    s->x = eig::ArrayXd::LinSpaced(n, 0.0, 1.0);
    s->r = s->x * s->R;
    s->rho              = 12.0 * (1.0 - 0.9 * s->x.square());
    s->P                = 3.0e17 * (1.0 - 0.95 * s->x.cube());
    s->Gamma_1          = 1.62 + 0.05 * (3.0 * s->x).sin();
    s->g                = 2.7e4 * s->x * (2.0 - s->x);
    s->dg_dr            = 2.7e4 * (2.0 - 2.0 * s->x) / s->R;
    s->dlnrho_dlnr      = -3.0 * s->x.square() / (1.0 - 0.9 * s->x.square() + 1e-3);
    s->dGamma1_dlnrho_s = 0.11 * (4.0 * s->x).cos();
    return s;
}

static Eigenfunction fake_ef(const std::shared_ptr<Star>& s, int l, int n_pg, double om) {
    Eigenfunction e;
    e.l = l;
    e.n_pg = n_pg;
    e.omega = om;
    const eig::ArrayXd& x = s->x;
    const double k = 3.0 + n_pg;
    e.xi_r          = 1.0e8 * (k * M_PI * x).sin() * (1.0 + x);
    e.xi_h          = 4.0e7 * (k * M_PI * x).cos() * x;
    e.div_xi        = 0.7 * (2.0 * k * x).sin() + 0.2;
    e.delta_Phi     = 1.0e14 * (k * x).cos();
    e.ddelta_Phi_dr = -1.0e14 * k * (k * x).sin() / s->R;
    e.starptr = s;
    return e;
}

int main() {
    /* --- angular --- */
    assert(satisfies_selection_rules(1, 2, 3));
    assert(!satisfies_selection_rules(1, 2, 2));          // odd parity
    assert(!satisfies_selection_rules(1, 2, 5));          // triangle
    assert(gaunt_T(1, 2, 2, 0, 0, 0) == 0.0);

    // GSL leaves ~5e-17 here where the exact 3j is zero; snapping keeps
    // m_combos from reporting a coupling that does not exist.
    assert(gaunt_T(2, 3, 3, 0, -2, 2) == 0.0);

    // T is invariant under simultaneous permutation of the (l, m) pairs.
    close(gaunt_T(1, 2, 3, 1, -1, 0), gaunt_T(2, 3, 1, -1, 0, 1), 1e-14, "T cyclic");
    close(gaunt_T(1, 2, 3, 1, -1, 0), gaunt_T(3, 2, 1, 0, -1, 1), 1e-14, "T swap");

    long n_m = 0;
    for (auto& L : l_multisets(4)) {
        assert(satisfies_selection_rules(L[0], L[1], L[2]));
        const auto ms = m_combos(L[0], L[1], L[2]);
        std::set<std::array<int, 3>> uniq(ms.begin(), ms.end());
        assert(uniq.size() == ms.size());                 // no combination emitted twice
        for (auto& m : ms) {
            assert(m[0] + m[1] + m[2] == 0);
            assert(gaunt_T(L[0], L[1], L[2], m[0], m[1], m[2]) != 0.0);
        }
        n_m += long(ms.size());
    }
    assert(n_m == 284);                                   // matches coupling/angular.py
    { // angular_factors on the trivial triplet
        AngularFactors f = angular_factors(0, 0, 0, 0, 0, 0);
        close(f.T, 1.0 / std::sqrt(4.0 * M_PI), 1e-14, "T(000)");
        for (double v : {f.F_a, f.F_b, f.F_c, f.G_a, f.G_b, f.G_c, f.S}) assert(v == 0.0);
    }

    /* --- kappa --- */
    auto star = fake_star(401);
    const Eigenfunction a = fake_ef(star, 2, 5, 1.4e-3);
    const Eigenfunction b = fake_ef(star, 2, 7, 1.9e-3);
    const Eigenfunction c = fake_ef(star, 2, 12, 3.1e-3);
    const EF3 efs = {&a, &b, &c};

    auto [bg, f] = on_grid(efs, 1);
    assert(bg.r.size() == star->x.size());

    // The one real restructure in the port: the seven basis densities
    // contracted with the angular factors must equal the eight group
    // densities summed, pointwise -- not merely in the integral.
    eig::ArrayXd (*const groups[8])(const Grid&, const F3&, const AngularFactors&) =
        {A55, A56, A57, A58, A59, A60, A61, A62};
    for (auto ms : {std::array<int,3>{0,0,0}, {1,-1,0}, {2,-1,-1}, {-2,1,1}}) {
        const AngularFactors ang = angular_factors(a.l, b.l, c.l, ms[0], ms[1], ms[2]);
        eig::ArrayXd from_groups = eig::ArrayXd::Zero(bg.r.size());
        for (int q = 0; q < 8; ++q) from_groups += groups[q](bg, f, ang);

        const auto d = basis_densities(bg, f);
        const auto af = ang.basis();
        eig::ArrayXd from_basis = eig::ArrayXd::Zero(bg.r.size());
        for (int j = 0; j < 7; ++j) from_basis += af[j] * d[j];

        const double scale = from_groups.abs().maxCoeff();
        assert(scale > 0.0);
        const double err = (from_basis - from_groups).abs().maxCoeff() / scale;
        if (err > 1e-12) {
            std::fprintf(stderr, "FAIL basis != groups for m=(%d,%d,%d): %.3e\n",
                         ms[0], ms[1], ms[2], err);
            assert(false);
        }
        // and the integral kappa_abc reports agrees with the group sum
        KappaResult k = kappa_abc(a, b, c, ms, 1);
        close(k.kappa, k.groups[0] + k.groups[1] + k.groups[2] + k.groups[3]
                     + k.groups[4] + k.groups[5] + k.groups[6] + k.groups[7],
              1e-10, "kappa == sum of groups");
    }

    // A56 splits exactly into its two pieces.
    { const AngularFactors ang = angular_factors(2, 2, 2, 0, 0, 0);
      auto [hor, rad] = A56_parts(bg, f, ang);
      assert((hor + rad - A56(bg, f, ang)).abs().maxCoeff()
             <= 1e-12 * A56(bg, f, ang).abs().maxCoeff()); }

    // Distinct Star objects holding identical grids take the union-merge branch.
    // Every node of the master grid is a node of all three, so the interpolants
    // reproduce them and the answer must match the shared-grid path.
    {
        Eigenfunction a2 = a, b2 = b, c2 = c;
        a2.starptr = std::make_shared<Star>(*star);
        b2.starptr = std::make_shared<Star>(*star);
        c2.starptr = std::make_shared<Star>(*star);
        assert(a2.starptr != b2.starptr);
        for (int refine : {1, 2, 4}) {
            KappaResult s = kappa_abc(a, b, c, {0, 0, 0}, refine);
            KappaResult u = kappa_abc(a2, b2, c2, {0, 0, 0}, refine);
            close(u.kappa, s.kappa, 1e-10, "union branch == shared branch");
        }
    }

    // Refinement settles, and the ladder stops once spline and Simpson agree.
    {
        KappaResult k8 = kappa_abc(a, b, c, {0, 0, 0}, 8);
        KappaResult lad = kappa_abc(a, b, c, {0, 0, 0});
        assert(lad.converged);
        assert(lad.quadrature_residual() <= QUAD_TOL);
        close(lad.kappa, k8.kappa, 1e-4, "ladder vs refine=8");
    }

    // kappa_all_m must reproduce kappa_abc mode by mode.
    {
        auto ms = m_combos(a.l, b.l, c.l);
        auto [ks, n] = kappa_all_m(a, b, c, ms, 2);
        assert(n == 2 && ks.size() == ms.size());
        for (size_t q = 0; q < ms.size(); ++q)
            close(ks[q], kappa_abc(a, b, c, ms[q], 2).kappa, 1e-11, "kappa_all_m vs kappa_abc");
    }

    std::printf("kappa: all checks passed\n");
    return 0;
}
