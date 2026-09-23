#include "angular.hpp"

#include <gsl/gsl_sf_coupling.h>

#include <cmath>


double wigner_3j (int l1, int l2, int l3, int m1, int m2, int m3) {
    return gsl_sf_coupling_3j(2*l1, 2*l2, 2*l3, 2*m1, 2*m2, 2*m3);
}

bool satisfies_selection_rules(int l_a, int l_b, int l_c) {
    return std::abs(l_b - l_c) <= l_a
            && l_a <= l_b + l_c
            && (l_a + l_b + l_c) % 2 == 0;
}

double gaunt_T (int l_a, int l_b, int l_c, int m_a, int m_b, int m_c) {
    if (m_a + m_b + m_c != 0 || !satisfies_selection_rules(l_a, l_b, l_c))
        return 0.0;

    const double norm = std::sqrt((2.0*l_a + 1) * (2.0*l_b + 1) * (2.0*l_c + 1)
                                 / (4.0 * M_PI));
    const double t = norm * wigner_3j(l_a, l_b, l_c, m_a, m_b, m_c)
                          * wigner_3j(l_a, l_b, l_c, 0, 0, 0);

    // GSL returns a small nonzero value instead of a hard zero
    return std::abs(t) < 1e-12 ? 0.0 : t;
}

AngularFactors angular_factors (int l_a, int l_b, int l_c, int m_a, int m_b, int m_c) {
    const double La = double(l_a) * (l_a + 1);
    const double Lb = double(l_b) * (l_b + 1);
    const double Lc = double(l_c) * (l_c + 1);

    AngularFactors a;
    a.T   = gaunt_T(l_a, l_b, l_c, m_a, m_b, m_c);
    a.F_a = 0.5 * a.T * (Lb + Lc - La);
    a.F_b = 0.5 * a.T * (Lc + La - Lb);
    a.F_c = 0.5 * a.T * (La + Lb - Lc);
    a.G_a = 0.25 * a.T * (La*La - (Lb - Lc)*(Lb - Lc));
    a.G_b = 0.25 * a.T * (Lb*Lb - (Lc - La)*(Lc - La));
    a.G_c = 0.25 * a.T * (Lc*Lc - (La - Lb)*(La - Lb));
    a.S   = 0.5 * (La * a.F_a + Lb * a.F_b + Lc * a.F_c);

    return a;
}

std::vector<std::array<int, 3>> l_multisets (int l_max) {
    std::vector<std::array<int, 3>> out;
    for (int a=0 ; a<=l_max ; a++)
        for (int b=a ; b<=l_max ; b++)
            for (int c=b ; c<=l_max ; c++)
                if (satisfies_selection_rules(a, b, c)) out.push_back({a, b, c});
    return out;
}

std::vector<std::array<int, 3>> m_combos (int l_a, int l_b, int l_c) {
    std::vector<std::array<int, 3>> out;
    for (int m_a=-l_a ; m_a<=l_a ; m_a++)
        for (int m_b=-l_b ; m_b<=l_b ; m_b++) {
            const int m_c = -m_a - m_b;
            if (std::abs(m_c) <= l_c && gaunt_T(l_a, l_b, l_c, m_a, m_b, m_c) != 0.0)
                out.push_back({m_a, m_b, m_c});
        }
    return out;
}
