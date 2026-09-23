#pragma once
#include <array>
#include <vector>

struct AngularFactors {
    double T, F_a, F_b, F_c, G_a, G_b, G_c, S;
    std::array<double, 7> basis() const {
        return {T, F_a, F_b, F_c, G_a, G_b, G_c}; 
    }
};

double wigner_3j(int l1, int l2, int l3, int m1, int m2, int m3);

bool satisfies_selection_rules(int l_a, int l_b, int l_c);

double gaunt_T(int l_a, int l_b, int l_c, int m_a, int m_b, int m_c);

AngularFactors angular_factors(int l_a, int l_b, int l_c, int m_a, int m_b, int m_c);

std::vector<std::array<int, 3>> l_multisets(int l_max);

std::vector<std::array<int, 3>> m_combos(int l_a, int l_b, int l_c);
