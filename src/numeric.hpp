#pragma once
#include <Eigen/Dense>

namespace eig = Eigen;

/* 
    Used manually written interpolation adn ingetration instead of GSL.
    GSL uses natural spline and steffen interpolation, whereas scipy uses
    not-a-knot and Fritsch-Carlson interpolation.
    Not exactly sure which one would be better, but kept the scipy default for now.
*/

namespace num {

struct CubicSpline {
    eig::ArrayXd x, y, M;                       // M = second derivative at each knot

    CubicSpline (const eig::ArrayXd& x, const eig::ArrayXd& y);

    double eval (double xq) const;
    eig::ArrayXd eval (const eig::ArrayXd& xq) const;
    eig::ArrayXd deriv_at_knots () const;
    eig::ArrayXd cumulative_integral () const;
};

struct Pchip {
    eig::ArrayXd x, y, d;                       // d = first derivative at each knot

    Pchip(const eig::ArrayXd& x, const eig::ArrayXd& y);

    double eval (double xq) const;
    eig::ArrayXd eval (const eig::ArrayXd& xq) const;
};

double trapz (const eig::ArrayXd& y, const eig::ArrayXd& x);
double simpson (const eig::ArrayXd& y, const eig::ArrayXd& x);

// xq need not be sorted; x must be ascending
eig::ArrayXd interp (const eig::ArrayXd& xq, const eig::ArrayXd& x, const eig::ArrayXd& y);

}
