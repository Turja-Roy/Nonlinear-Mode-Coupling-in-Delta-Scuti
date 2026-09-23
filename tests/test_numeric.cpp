// g++ -std=c++17 -I/usr/include/eigen3 -Isrc tests/test_numeric.cpp src/numeric.cpp -o t && ./t
#include "numeric.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>

static void close(double a, double b, double tol, const char* what) {
    if (std::abs(a - b) > tol * (1.0 + std::abs(b))) {
        std::fprintf(stderr, "FAIL %s: %.17g != %.17g\n", what, a, b);
        assert(false);
    }
}

static eig::ArrayXd arr(std::initializer_list<double> v) {
    eig::ArrayXd a(v.size());
    int i = 0;
    for (double d : v) a[i++] = d;
    return a;
}

int main() {
    // Not-a-knot is exact for cubics: value, slope and integral all reproduce.
    const eig::ArrayXd x = arr({0.0, 0.3, 0.7, 1.4, 2.0, 3.1, 3.5, 5.0});
    auto p  = [](double t) { return 2*t*t*t - 3*t*t + t + 5; };
    auto dp = [](double t) { return 6*t*t - 6*t + 1; };
    auto ip = [](double t) { return 0.5*t*t*t*t - t*t*t + 0.5*t*t + 5*t; };

    eig::ArrayXd y = x.unaryExpr(p);
    num::CubicSpline s(x, y);

    for (double q : {0.15, 0.7, 1.0, 2.5, 4.9, -0.4, 5.6})       // last two extrapolate
        close(s.eval(q), p(q), 1e-11, "spline eval");

    eig::ArrayXd ds = s.deriv_at_knots(), cum = s.cumulative_integral();
    for (eig::Index i = 0; i < x.size(); ++i) {
        close(ds[i], dp(x[i]), 1e-11, "spline deriv");
        close(cum[i], ip(x[i]) - ip(x[0]), 1e-11, "spline cumulative");
    }
    assert(cum[0] == 0.0);

    // Uniform grid: h0 == h1 exactly, which is what subdivide() produces and
    // what degenerates the naive not-a-knot elimination into a zero pivot.
    for (int n : {4, 5, 9, 64}) {
        eig::ArrayXd xu(n);
        for (int i = 0; i < n; ++i) xu[i] = 0.25 * i;
        num::CubicSpline su(xu, xu.unaryExpr(p));
        assert(std::isfinite(su.M.sum()));
        close(su.eval(0.3), p(0.3), 1e-11, "uniform spline eval");
        close(su.deriv_at_knots()[0], dp(0.0), 1e-11, "uniform spline deriv at x0");
        close(su.cumulative_integral()[n-1], ip(xu[n-1]) - ip(0.0), 1e-11, "uniform cumulative");
    }

    // Degenerate sizes: 3 points collapse to the parabola, 2 to a line.
    const eig::ArrayXd x3 = arr({0.0, 1.0, 3.0});
    num::CubicSpline s3(x3, x3.unaryExpr([](double t) { return t*t - 2*t + 1; }));
    close(s3.eval(2.0), 1.0, 1e-12, "n=3 parabola");
    close(s3.eval(-1.0), 4.0, 1e-12, "n=3 parabola extrapolated");

    num::CubicSpline s2(arr({1.0, 4.0}), arr({2.0, 8.0}));
    close(s2.eval(2.0), 4.0, 1e-12, "n=2 line");
    close(s2.deriv_at_knots()[0], 2.0, 1e-12, "n=2 slope");

    // PCHIP: exact on a line, monotone-preserving, flat at an interior extremum.
    const eig::ArrayXd xp = arr({0.0, 1.0, 2.5, 3.0, 5.0});
    num::Pchip lin(xp, 3.0 * xp + 1.0);
    for (double q : {0.5, 2.0, 4.4}) close(lin.eval(q), 3*q + 1, 1e-12, "pchip linear");

    const eig::ArrayXd ym = arr({0.0, 0.0, 0.0, 1.0, 1.0});      // classic overshoot trap
    num::Pchip mono(xp, ym);
    double prev = -1e300;
    for (int k = 0; k <= 200; ++k) {
        double v = mono.eval(5.0 * k / 200.0);
        assert(v >= prev - 1e-12 && v >= -1e-12 && v <= 1.0 + 1e-12);
        prev = v;
    }

    num::Pchip peak(xp, arr({0.0, 1.0, 2.0, 1.0, 0.0}));
    close(peak.d[2], 0.0, 1e-12, "pchip flat at extremum");

    // trapz exact on a line; non-uniform Simpson exact on a quadratic.
    close(num::trapz(3.0 * x + 1.0, x), 1.5*25 + 5, 1e-11, "trapz linear");
    auto q2 = [](double t) { return 4*t*t - t + 2; };
    {
        eig::ArrayXd xe = arr({0.0, 0.4, 1.0, 1.1, 2.0});              // 4 cells: exact
        close(num::simpson(xe.unaryExpr(q2), xe), (4.0/3)*8 - 2 + 4, 1e-11, "simpson even");

        eig::ArrayXd xo = arr({0.0, 0.4, 1.0, 1.1, 2.0, 3.0});         // 5 cells: tail is
        double exact = (4.0/3)*27 - 4.5 + 6;                           // one trapezoid, so
        double h = 1.0, err = h*h*h * 8.0 / 12.0;                      // it overshoots by
        close(num::simpson(xo.unaryExpr(q2), xo), exact + err, 1e-11, "simpson odd tail");
    }

    // np.interp clamps rather than extrapolating -- load-bearing in mesa_on_gyre_grid.
    eig::ArrayXd g = num::interp(arr({-5.0, 0.5, 2.0, 99.0}), arr({0.0, 1.0, 3.0}),
                                                              arr({10.0, 20.0, 40.0}));
    close(g[0], 10.0, 1e-12, "interp clamp low");
    close(g[1], 15.0, 1e-12, "interp mid");
    close(g[2], 30.0, 1e-12, "interp mid");
    close(g[3], 40.0, 1e-12, "interp clamp high");

    std::printf("numeric: all checks passed\n");
    return 0;
}
