#include "numeric.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace num {

namespace {

// Last cell whose left knot is <= xq, clamped so out-of-range queries
// extrapolate with the end cell's polynomial (scipy's extrapolate=True).
eig::Index cell (const eig::ArrayXd& x, double xq) {
    const double* p = std::upper_bound(x.data(), x.data() + x.size(), xq);
    return std::clamp<eig::Index>( 
        (p-x.data())-1 , 0 , x.size()-2
    );
}

int sgn (double v) {
    return (v > 0) - (v < 0); 
}

}


/* ------------------------
    Not-a-knot cubic spline
   ------------------------
*/

CubicSpline::CubicSpline (const eig::ArrayXd& x_, const eig::ArrayXd& y_) : x(x_), y(y_) {
    const eig::Index n = x.size();

    if (n != y.size())
        throw std::invalid_argument("CubicSpline: size mismatch");
    if (n < 2)
        throw std::invalid_argument("CubicSpline: needs >= 2 points");

    M = eig::ArrayXd::Zero(n);
    if (n == 2) return;             // straight line, M = 0

    const eig::ArrayXd h = x.tail(n-1) - x.head(n-1);           // n-1 cells
    const eig::ArrayXd m = (y.tail(n-1) - y.head(n-1)) / h;     // secants

    if (n == 3) {
        M.setConstant( 2.0 * (m[1]-m[0]) / (x[2]-x[0]) );   // collapse to the
        return;                                             // parabola through the 3
    }

    /* Unknowns are M[1..n-2]. Each not-a-knot row gives one end knot in terms
       of its two neighbours; substituting those into the first and last
       interior rows leaves a tridiagonal, strictly diagonally dominant system.
       Eliminating the other way -- making the not-a-knot rows themselves
       tridiagonal -- puts (h1^2 - h0^2)/h1 on the pivot, which is exactly zero
       on a uniform grid, and subdivide() produces uniform cells. */
    const eig::Index K = n - 2;
    eig::ArrayXd a(K), b(K), c(K), d(K);
    for (eig::Index t=0 ; t<K ; t++) {
        const eig::Index i = t + 1;
        a[t] = h[i-1];
        b[t] = 2.0 * (h[i-1] + h[i]);
        c[t] = h[i];
        d[t] = 6.0 * (m[i] - m[i-1]);
    }

    a[0] = 0.0;                                     // M0 = ((h0+h1) M1 - h0 M2) / h1
    b[0] = (h[0] + h[1]) * (h[0] + 2.0*h[1]) / h[1];
    c[0] = (h[1]*h[1] - h[0]*h[0]) / h[1];

    const double p = h[n-3], q = h[n-2];            // M_{n-1}, mirrored
    a[K-1] = (p*p - q*q) / p;
    b[K-1] = (p + q) * (2.0*p + q) / p;
    c[K-1] = 0.0;

    c[0] /= b[0];
    d[0] /= b[0];
    for (eig::Index t=1 ; t<K ; t++) {
        const double den = b[t] - a[t] * c[t-1];
        c[t] /= den;
        d[t] = (d[t] - a[t] * d[t-1]) / den;
    }
    M[n-2] = d[K-1];
    for (eig::Index t=K-2 ; t>=0 ; t--)
        M[t+1] = d[t] - c[t] * M[t+2];

    M[0]   = ((h[0] + h[1]) * M[1] - h[0] * M[2]) / h[1];
    M[n-1] = ((p + q) * M[n-2] - q * M[n-3]) / p;
}

double CubicSpline::eval (double xq) const {
    const eig::Index i = cell(x, xq);
    const double h = x[i + 1] - x[i], t = (xq - x[i]) / h;
    return y[i] * (1.0 - t) + y[i + 1] * t
         - (h * h / 6.0) * t * (1.0 - t) * ((2.0 - t) * M[i] + (1.0 + t) * M[i + 1]);
}

eig::ArrayXd CubicSpline::eval (const eig::ArrayXd& xq) const {
    eig::ArrayXd out(xq.size());
    for (eig::Index k = 0; k < xq.size(); ++k) out[k] = eval(xq[k]);
    return out;
}

eig::ArrayXd CubicSpline::deriv_at_knots() const {
    const eig::Index n = x.size();
    eig::ArrayXd out(n);
    for (eig::Index i = 0; i < n - 1; ++i) {
        const double h = x[i + 1] - x[i];
        out[i] = (y[i + 1] - y[i]) / h - h * (2.0 * M[i] + M[i + 1]) / 6.0;
    }
    const double h = x[n - 1] - x[n - 2];                 // right end of the last cell
    out[n - 1] = (y[n - 1] - y[n - 2]) / h + h * (M[n - 2] + 2.0 * M[n - 1]) / 6.0;
    return out;
}

eig::ArrayXd CubicSpline::cumulative_integral() const {
    const eig::Index n = x.size();
    eig::ArrayXd out(n);
    out[0] = 0.0;
    for (eig::Index i = 0; i < n - 1; ++i) {
        const double h = x[i + 1] - x[i];
        out[i + 1] = out[i] + h * (y[i] + y[i + 1]) / 2.0
                            - h * h * h * (M[i] + M[i + 1]) / 24.0;
    }
    return out;
}


/* ----------------------
    Fritsch-Carlson PCHIP
   ----------------------
*/

Pchip::Pchip (const eig::ArrayXd& x_, const eig::ArrayXd& y_) : x(x_), y(y_) {
    const eig::Index n = x.size();
    if (n != y.size())
        throw std::invalid_argument("Pchip: size mismatch");
    if (n < 2)
        throw std::invalid_argument("Pchip: needs >= 2 points");

    const eig::ArrayXd h = x.tail(n-1) - x.head(n-1);
    const eig::ArrayXd m = (y.tail(n-1) - y.head(n-1)) / h;

    d = eig::ArrayXd::Zero(n);
    if (n == 2) { d.setConstant(m[0]); return; }

    for (eig::Index i=1 ; i<n-1 ; i++) {
        if (m[i-1] * m[i] <= 0.0) continue;             // extremum: flat, no overshoot
        const double w1 = 2.0 * h[i] + h[i-1];
        const double w2 = h[i] + 2.0 * h[i-1];
        d[i] = (w1+w2) / (w1 / m[i-1] + w2 / m[i]);     // weighted harmonic mean
    }

    // scipy's _edge_case: one-sided quadratic slope, then two shape limiters.
    auto edge = [](double h0, double h1, double m0, double m1) {
        double v = ( (2.0*h0 + h1)*m0 - h0*m1 ) / (h0+h1);
        if (sgn(v) != sgn(m0))
            return 0.0;
        if (sgn(m0) != sgn(m1) && std::abs(v) > 3.0 * std::abs(m0))
            return 3.0 * m0;
        return v;
    };
    d[0] = edge(h[0], h[1], m[0], m[1]);
    d[n-1] = edge(h[n-2], h[n-3], m[n-2], m[n-3]);
}

double Pchip::eval (double xq) const {
    const eig::Index i = cell(x, xq);
    const double h = x[i + 1] - x[i], t = (xq - x[i]) / h;
    const double t2 = t * t, t3 = t2 * t;

    return (2.0 * t3 - 3.0 * t2 + 1.0) * y[i]                     // cubic Hermite
         + (t3 - 2.0 * t2 + t) * h * d[i]
         + (-2.0 * t3 + 3.0 * t2) * y[i + 1]
         + (t3 - t2) * h * d[i + 1];
}

eig::ArrayXd Pchip::eval (const eig::ArrayXd& xq) const {
    eig::ArrayXd out(xq.size());
    for (eig::Index k=0 ; k<xq.size() ; k++) out[k] = eval(xq[k]);
    return out;
}


/* -----------------------
    Quadrature and interp
   -----------------------
*/

double trapz (const eig::ArrayXd& y, const eig::ArrayXd& x) {
    const eig::Index n = x.size();
    if (n != y.size())
        throw std::invalid_argument("trapz: size mismatch");
    if (n < 2)
        return 0.0;
    return 0.5 * ( (x.tail(n-1) - x.head(n-1)) * (y.tail(n-1) + y.head(n-1)) ).sum();
}

double simpson(const eig::ArrayXd& y, const eig::ArrayXd& x) {
    const eig::Index n = x.size();
    if (n != y.size()) throw std::invalid_argument("simpson: size mismatch");
    if (n < 3) return trapz(y, x);

    double s = 0.0;
    eig::Index i = 0;
    for (; i + 2 < n; i += 2) {
        const double h0 = x[i + 1] - x[i], h1 = x[i + 2] - x[i + 1];
        s += (h0 + h1) / 6.0 * ((2.0 - h1 / h0) * y[i]
                              + (h0 + h1) * (h0 + h1) / (h0 * h1) * y[i + 1]
                              + (2.0 - h0 / h1) * y[i + 2]);
    }

    if (i+1 < n) s += 0.5 * (x[i+1] - x[i]) * (y[i] + y[i+1]);
    return s;
}

eig::ArrayXd interp (const eig::ArrayXd& xq, const eig::ArrayXd& x, const eig::ArrayXd& y) {
    const eig::Index n = x.size();
    if (n != y.size())
        throw std::invalid_argument("interp: size mismatch");
    if (n == 0)
        throw std::invalid_argument("interp: empty table");

    eig::ArrayXd out(xq.size());
    for (eig::Index k=0 ; k<xq.size() ; k++) {
        const double q = xq[k];
        if (n == 1 || q <= x[0]) {
            out[k] = y[0];
            continue; 
        }
        if (q >= x[n - 1]) {
            out[k] = y[n - 1];
            continue; 
        }
        const eig::Index i = cell(x, q);
        out[k] = y[i] + (y[i+1] - y[i]) * (q - x[i]) / (x[i+1] - x[i]);
    }
    return out;
}

}
