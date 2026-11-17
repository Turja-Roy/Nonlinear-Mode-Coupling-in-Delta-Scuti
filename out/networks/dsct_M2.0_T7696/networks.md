# dsct_M2.0_T7696: ladders and runs

Written by scripts/network_index.py.

Notation. (l,n): degree and radial order of a mode; f: frequency; gamma: linear rate, < 0 driven (parent),
> 0 damped (daughter). N: number of daughters. Pair p -> b + c: parametric triplet; harmonic x + x -> c and
sum x + y -> c: direct triplets. Delta: detuning; kappa: coupling coefficient. E_th: parent energy above
which a pair grows (Weinberg+12 eq. 50); E_eq: parent energy at a triplet's fixed point (MW25 eq. 6), for a
sum triplet E_x + E_y; for pairs recomputed here from the network file. gamma sum: gamma_p + gamma_b + gamma_c
of a pair, gamma_x + gamma_y + gamma_c of a direct triplet (y = x for a harmonic); <= 0: a pair's fixed point
is unstable, a harmonic triplet cannot hold x alone (a sum triplet's fixed point is unstable at any value).
Energies in E_star = GM^2/R. Shared daughter: coupled to both parents. Run columns,
second half of the run: <E_par>, mean summed parent energy; swing, (max - min)/mean of it; n_active,
daughters above 1e-3 of the largest. <D>_T/<P>_T: damping power D = sum 2 gamma_d E_d over driving power
P = sum 2|gamma_p| E_p, each a moving time average (1/T) * integral over (t - T, t], window T = t_stop/10,
taken at t = t_stop. stop: end (t = 1e13 s), runaway (summed energy 1e-2), wall (time limit).

## ladder_1p

### ladder_1p/l0_n+1

Parents: (0,+1) f 10.715 c/d, gamma -7.904e-11 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (0,+1) -> (6,-55) + (6,-10) | +3.995e-03 | +1.288e-02 | 4.814e-07 | 4.971e-07 | +4.75e-09 |  |
| pair | (0,+1) -> (4,-32) + (4,-7) | +2.463e-02 | +1.655e-01 | 1.465e-06 | 8.393e-06 | +5.67e-11 |  |
| pair | (0,+1) -> (5,-41) + (5,-9) | -3.058e-02 | +1.201e-01 | 1.550e-06 | 2.021e-06 | +5.55e-10 |  |
| pair | (0,+1) -> (6,-38) + (6,-11) | +1.378e-02 | +8.414e-02 | 1.609e-06 | 2.465e-06 | +3.32e-10 |  |
| pair | (0,+1) -> (6,-56) + (6,-10) | -2.655e-02 | +4.178e-02 | 2.520e-06 | 2.626e-06 | +3.85e-09 |  |
| pair | (0,+1) -> (4,-33) + (4,-7) | -4.778e-02 | +2.058e-01 | 2.595e-06 | 7.431e-06 | +1.14e-10 |  |
| pair | (0,+1) -> (3,-20) + (3,-6) | -2.760e-02 | +3.699e-01 | 3.064e-06 | 8.638e-08 | -6.77e-11 |  |
| pair | (0,+1) -> (4,-23) + (4,-8) | -1.395e-02 | +1.398e-01 | 4.643e-06 | 5.085e-07 | -5.94e-11 |  |
| pair | (0,+1) -> (6,-54) + (6,-10) | +4.256e-02 | +4.230e-02 | 5.313e-06 | 5.506e-06 | +4.43e-09 |  |
| pair | (0,+1) -> (2,-15) + (2,-4) | -2.126e-02 | +2.692e-01 | 5.931e-06 | 1.608e-08 | -7.51e-11 |  |
| pair | (0,+1) -> (5,-39) + (5,-9) | +6.029e-02 | +1.307e-01 | 7.966e-06 | 1.258e-05 | +3.09e-10 |  |
| pair | (0,+1) -> (5,-40) + (5,-9) | +2.061e-02 | +3.151e-02 | 1.534e-05 | 2.348e-05 | +3.33e-10 |  |
| pair | (0,+1) -> (5,-23) + (5,-11) | -3.075e-02 | +2.191e-01 | 1.615e-05 | 2.397e-06 | -5.71e-11 |  |
| pair | (0,+1) -> (6,-39) + (6,-11) | -5.254e-02 | +1.085e-01 | 1.655e-05 | 2.729e-05 | +2.79e-10 |  |
| pair | (0,+1) -> (4,-31) + (4,-7) | +7.947e-02 | +1.929e-01 | 1.888e-05 | 3.510e-01 | -5.76e-13 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | end | 1e+13 | 5.0899e-07 | 184 | 0.7497 | 1 |
| d6 | 6 | 3 |  |  | end | 1e+13 | 1.0925e-06 | 2.08 | 1.122 | 3 |
| d8 | 8 | 4 |  |  | end | 1e+13 | 1.0927e-06 | 2.08 | 1.122 | 3 |
| d9 | 9 | 5 |  |  | end | 1e+13 | 4.7526e-07 | 4.45 | 1.071 | 3 |
| d10 | 10 | 6 |  |  | end | 1e+13 | 3.9981e-07 | 0.172 | 1.001 | 1 |
| d12 | 12 | 7 |  |  | end | 1e+13 | 4.0724e-07 | 31.2 | 0.6744 | 1 |
| d14 | 14 | 8 |  |  | end | 1e+13 | 4.0724e-07 | 31.2 | 0.6744 | 1 |
| d15 | 15 | 9 |  |  | end | 1e+13 | 3.8832e-07 | 22.4 | 0.9654 | 1 |
| d17 | 17 | 10 |  |  | end | 1e+13 | 3.8836e-07 | 22.3 | 0.9651 | 1 |
| d18 | 18 | 11 |  |  | end | 1e+13 | 3.8836e-07 | 22.3 | 0.9651 | 1 |
| d19 | 19 | 12 |  |  | end | 1e+13 | 3.8836e-07 | 22.3 | 0.9651 | 1 |
| d21 | 21 | 13 |  |  | end | 1e+13 | 3.8836e-07 | 22.3 | 0.9651 | 1 |
| d22 | 22 | 14 |  |  | end | 1e+13 | 3.8836e-07 | 22.3 | 0.9651 | 1 |
| d23 | 23 | 15 |  |  | end | 1e+13 | 3.8836e-07 | 22.3 | 0.9651 | 1 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 17 on.

### ladder_1p/l0_n+1_ladder2

Parents: (0,+1) f 10.715 c/d, gamma -7.904e-11 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first; only pairs with an E_eq.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (0,+1) -> (6,-55) + (6,-10) | +3.995e-03 | +1.288e-02 | 4.814e-07 | 4.971e-07 | +4.75e-09 |  |
| pair | (0,+1) -> (4,-32) + (4,-7) | +2.463e-02 | +1.655e-01 | 1.465e-06 | 8.393e-06 | +5.67e-11 |  |
| pair | (0,+1) -> (5,-41) + (5,-9) | -3.058e-02 | +1.201e-01 | 1.550e-06 | 2.021e-06 | +5.55e-10 |  |
| pair | (0,+1) -> (6,-38) + (6,-11) | +1.378e-02 | +8.414e-02 | 1.609e-06 | 2.465e-06 | +3.32e-10 |  |
| pair | (0,+1) -> (6,-56) + (6,-10) | -2.655e-02 | +4.178e-02 | 2.520e-06 | 2.626e-06 | +3.85e-09 |  |
| pair | (0,+1) -> (4,-33) + (4,-7) | -4.778e-02 | +2.058e-01 | 2.595e-06 | 7.431e-06 | +1.14e-10 |  |
| pair | (0,+1) -> (6,-54) + (6,-10) | +4.256e-02 | +4.230e-02 | 5.313e-06 | 5.506e-06 | +4.43e-09 |  |
| pair | (0,+1) -> (5,-39) + (5,-9) | +6.029e-02 | +1.307e-01 | 7.966e-06 | 1.258e-05 | +3.09e-10 |  |
| pair | (0,+1) -> (5,-40) + (5,-9) | +2.061e-02 | +3.151e-02 | 1.534e-05 | 2.348e-05 | +3.33e-10 |  |
| pair | (0,+1) -> (6,-39) + (6,-11) | -5.254e-02 | +1.085e-01 | 1.655e-05 | 2.729e-05 | +2.79e-10 |  |
| pair | (0,+1) -> (5,-29) + (5,-10) | -5.395e-02 | +1.706e-01 | 2.016e-05 | 5.182e-03 | +5.26e-12 |  |
| pair | (0,+1) -> (6,-58) + (6,-10) | -8.595e-02 | +3.329e-02 | 2.368e-05 | 2.422e-05 | +7.06e-09 |  |
| pair | (0,+1) -> (6,-32) + (6,-12) | -5.173e-02 | +1.452e-01 | 2.380e-05 | 1.698e-04 | +4.73e-11 |  |
| pair | (0,+1) -> (4,-34) + (4,-7) | -1.177e-01 | +1.515e-01 | 2.529e-05 | 5.884e-05 | +1.50e-10 |  |
| pair | (0,+1) -> (5,-43) + (5,-9) | -1.357e-01 | +1.171e-01 | 3.412e-05 | 4.467e-05 | +5.48e-10 |  |
| pair | (0,+1) -> (6,-37) + (6,-11) | +9.385e-02 | +1.215e-01 | 4.144e-05 | 6.986e-05 | +2.65e-10 |  |
| pair | (0,+1) -> (4,-30) + (4,-7) | +1.400e-01 | +1.939e-01 | 4.433e-05 | 9.930e-04 | +2.12e-11 |  |
| pair | (0,+1) -> (6,-52) + (6,-10) | +1.136e-01 | +5.291e-02 | 4.721e-05 | 5.077e-05 | +2.15e-09 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d23 | 23 | 17 |  |  | end | 1e+13 | 3.7241e-07 | 0.78 | 0.9793 | 1 |
| d24 | 24 | 18 |  |  | end | 1e+13 | 3.7325e-07 | 0.853 | 0.9846 | 1 |

Same run from some N on: not reached.

### ladder_1p/l0_n+2

Parents: (0,+2) f 13.832 c/d, gamma -1.906e-09 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (0,+2) -> (6,-24) + (6,-9) | +7.283e-02 | +9.585e-01 | 7.165e-07 | 1.858e-10 | -1.88e-09 |  |
| pair | (0,+2) -> (6,-25) + (6,-9) | -1.030e-01 | +9.688e-01 | 1.200e-06 | 4.709e-10 | -1.87e-09 |  |
| pair | (0,+2) -> (5,-18) + (5,-8) | -3.607e-02 | +2.975e-01 | 8.026e-06 | 1.133e-10 | -1.90e-09 |  |
| pair | (0,+2) -> (6,-26) + (6,-9) | -2.422e-01 | +8.753e-01 | 9.301e-06 | 2.971e-09 | -1.87e-09 |  |
| pair | (0,+2) -> (6,-23) + (6,-9) | +2.615e-01 | +8.251e-01 | 1.606e-05 | 2.241e-09 | -1.88e-09 |  |
| pair | (0,+2) -> (6,-27) + (6,-9) | -3.363e-01 | +8.553e-01 | 1.643e-05 | 7.293e-09 | -1.87e-09 |  |
| pair | (0,+2) -> (6,-14) + (6,-13) | +2.594e-02 | +2.099e-01 | 1.918e-05 | 2.338e-10 | -1.90e-09 |  |
| pair | (0,+2) -> (5,-13) + (5,-10) | +2.635e-02 | +1.315e-01 | 5.045e-05 | 2.255e-10 | -1.90e-09 |  |
| pair | (0,+2) -> (4,-12) + (4,-7) | +5.479e-02 | +2.230e-01 | 5.796e-05 | 1.447e-10 | -1.90e-09 |  |
| pair | (0,+2) -> (6,-19) + (6,-10) | +1.171e-01 | +3.386e-01 | 7.580e-05 | 3.166e-09 | -1.89e-09 |  |
| pair | (0,+2) -> (5,-17) + (5,-8) | +1.252e-01 | +3.444e-01 | 7.796e-05 | 8.521e-10 | -1.90e-09 |  |
| pair | (0,+2) -> (6,-20) + (6,-10) | -1.529e-01 | +3.816e-01 | 9.092e-05 | 5.605e-09 | -1.89e-09 |  |
| pair | (0,+2) -> (6,-22) + (6,-9) | +4.172e-01 | +6.291e-01 | 1.111e-04 | 5.518e-09 | -1.89e-09 |  |
| pair | (0,+2) -> (5,-19) + (5,-8) | -2.652e-01 | +4.626e-01 | 1.363e-04 | 3.968e-09 | -1.90e-09 |  |
| pair | (0,+2) -> (5,-15) + (5,-9) | +1.093e-01 | +2.378e-01 | 1.759e-04 | 1.674e-09 | -1.90e-09 |  |
| pair | (0,+2) -> (6,-16) + (6,-11) | +7.088e-02 | +1.599e-01 | 1.971e-04 | 3.927e-09 | -1.90e-09 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | runaway | 8.35e+09 | 6.3788e-04 | 8.45 | 0.007926 | 2 |
| d3 | 3 | 2 |  |  | runaway | 9.15e+09 | 5.7974e-04 | 13 | 0.007547 | 3 |
| d5 | 5 | 3 |  |  | runaway | 1.06e+10 | 5.0406e-04 | 10.8 | 0.009715 | 5 |
| d6 | 6 | 4 |  |  | runaway | 9.48e+09 | 5.6426e-04 | 8.45 | 0.009787 | 6 |
| d7 | 7 | 5 |  |  | runaway | 1.08e+10 | 5.0487e-04 | 13.7 | 0.009039 | 7 |
| d8 | 8 | 6 |  |  | runaway | 1.08e+10 | 3.6287e-04 | 8.18 | 0.01053 | 8 |
| d10 | 10 | 7 |  |  | runaway | 1.11e+10 | 2.1166e-04 | 8.3 | 0.01087 | 10 |
| d12 | 12 | 8 |  |  | runaway | 1.14e+10 | 4.5039e-04 | 15 | 0.008098 | 12 |
| d14 | 14 | 9 |  |  | runaway | 1.16e+10 | 4.5872e-04 | 13.3 | 0.007967 | 14 |
| d16 | 16 | 10 |  |  | runaway | 1.31e+10 | 4.0925e-04 | 13.4 | 0.009042 | 16 |
| d17 | 17 | 11 |  |  | runaway | 1.33e+10 | 3.8778e-04 | 13.1 | 0.009964 | 17 |
| d18 | 18 | 12 |  |  | runaway | 1.41e+10 | 3.8093e-04 | 11.7 | 0.01006 | 18 |
| d19 | 19 | 13 |  |  | runaway | 1.64e+10 | 3.2223e-04 | 15.3 | 0.00933 | 19 |
| d20 | 20 | 14 |  |  | runaway | 1.62e+10 | 3.2129e-04 | 11.5 | 0.01027 | 20 |
| d22 | 22 | 15 |  |  | runaway | 1.73e+10 | 2.9206e-04 | 19.5 | 0.009435 | 22 |
| d24 | 24 | 16 |  |  | runaway | 1.78e+10 | 2.6182e-04 | 15.3 | 0.01024 | 24 |

Same run from some N on: not reached.

### ladder_1p/l0_n+3

Parents: (0,+3) f 17.180 c/d, gamma -1.432e-08 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (0,+3) -> (6,-12) + (6,-9) | +6.164e-02 | +9.035e-02 | 2.664e-04 | 1.106e-11 | -1.43e-08 |  |
| pair | (0,+3) -> (5,-9) + (5,-9) | -8.472e-02 | -8.939e-02 | 7.683e-04 | 1.846e-11 | -1.43e-08 |  |
| pair | (0,+3) -> (4,-7) + (4,-7) | -1.414e-01 | -1.408e-01 | 8.686e-04 | 7.907e-12 | -1.43e-08 |  |
| pair | (0,+3) -> (6,-13) + (6,-9) | -2.591e-01 | +1.086e-01 | 3.096e-03 | 1.655e-10 | -1.43e-08 |  |
| pair | (0,+3) -> (6,-11) + (6,-10) | -3.203e-01 | +4.483e-02 | 4.410e-02 | 3.047e-09 | -1.43e-08 |  |
| pair | (0,+3) -> (5,-10) + (5,-8) | -2.949e-01 | -1.412e-02 | 3.348e-01 | 8.727e-09 | -1.43e-08 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | runaway | 8.82e+08 | 7.9220e-04 | 7.39 | 0.0001079 | 2 |
| d3 | 3 | 2 |  |  | runaway | 1.02e+09 | 6.8992e-04 | 11.2 | 0.0002437 | 3 |
| d4 | 4 | 3 |  |  | runaway | 9.74e+08 | 7.1314e-04 | 11.5 | 0.0001508 | 4 |
| d5 | 5 | 4 |  |  | runaway | 9.41e+08 | 7.4572e-04 | 8.69 | 0.0002027 | 5 |
| d7 | 7 | 5 |  |  | runaway | 9.41e+08 | 7.4848e-04 | 6.67 | 0.0002037 | 5 |
| d9 | 9 | 6 |  |  | runaway | 9.39e+08 | 7.4375e-04 | 8.07 | 0.0001985 | 5 |

Same run from some N on: not reached.

### ladder_1p/l0_n+4

Parents: (0,+4) f 20.577 c/d, gamma -6.490e-08 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| harmonic | (0,+4) + (0,+4) -> (0,+10) | +2.199e-01 | +1.431e+01 |  | 3.186e-09 | +3.05e-06 | 2.640e+03 |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  | (0,+4) + (0,+4) -> (0,+10) | end | 1e+13 | 3.1859e-09 | 2.13e-10 | 1.005 | 1 |

Same run from some N on: not reached.

### ladder_1p/l0_n+8

Parents: (0,+8) f 34.144 c/d, gamma -4.485e-07 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| harmonic | (0,+8) + (0,+8) -> (0,+17) | -2.924e-01 | +1.196e+02 |  | 5.758e-11 | +2.68e-05 | 1.695e+04 |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  | (0,+8) + (0,+8) -> (0,+17) | end | 1e+13 | 5.7585e-11 | 9.11e-11 | 0.9957 | 1 |

Same run from some N on: not reached.

### ladder_1p/l1_n+1

Parents: (1,+1) f 14.161 c/d, gamma -2.355e-09 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (1,+1) -> (6,-9) + (7,-26) | -8.827e-03 | +6.457e-01 | 2.096e-08 | 3.869e-12 | -2.32e-09 |  |
| pair | (1,+1) -> (5,-19) + (6,-9) | -4.469e-03 | +5.008e-01 | 2.667e-08 | 4.959e-13 | -2.35e-09 |  |
| pair | (1,+1) -> (4,-9) + (5,-11) | +8.258e-03 | +2.029e-01 | 1.865e-06 | 3.311e-12 | -2.35e-09 |  |
| pair | (1,+1) -> (6,-9) + (7,-27) | -1.179e-01 | +7.438e-01 | 2.137e-06 | 7.362e-10 | -2.31e-09 |  |
| pair | (1,+1) -> (6,-9) + (7,-25) | +1.381e-01 | +7.587e-01 | 3.019e-06 | 7.997e-10 | -2.32e-09 |  |
| pair | (1,+1) -> (5,-11) + (6,-13) | -1.294e-02 | +2.545e-01 | 3.169e-06 | 1.292e-11 | -2.35e-09 |  |
| pair | (1,+1) -> (5,-9) + (6,-17) | -1.599e-02 | +1.993e-01 | 4.666e-06 | 3.762e-11 | -2.35e-09 |  |
| pair | (1,+1) -> (6,-9) + (7,-28) | -2.672e-01 | +8.630e-01 | 5.749e-06 | 4.396e-09 | -2.29e-09 |  |
| pair | (1,+1) -> (6,-9) + (7,-29) | -4.203e-01 | +8.936e-01 | 1.120e-05 | 1.317e-08 | -2.28e-09 |  |
| pair | (1,+1) -> (6,-9) + (7,-24) | +3.334e-01 | +7.279e-01 | 2.158e-05 | 4.079e-09 | -2.32e-09 |  |
| pair | (1,+1) -> (5,-8) + (6,-20) | -1.294e-01 | +3.676e-01 | 3.470e-05 | 1.237e-09 | -2.34e-09 |  |
| pair | (1,+1) -> (4,-10) + (5,-10) | +3.120e-02 | +1.836e-01 | 3.608e-05 | 6.049e-11 | -2.35e-09 |  |
| pair | (1,+1) -> (5,-20) + (6,-9) | -2.406e-01 | +6.136e-01 | 4.369e-05 | 1.299e-09 | -2.34e-09 |  |
| pair | (1,+1) -> (5,-8) + (6,-19) | +1.407e-01 | +3.175e-01 | 6.329e-05 | 1.487e-09 | -2.34e-09 |  |
| pair | (1,+1) -> (4,-8) + (5,-12) | -5.099e-02 | +1.669e-01 | 9.380e-05 | 1.819e-10 | -2.35e-09 |  |
| pair | (1,+1) -> (5,-15) + (6,-10) | +5.586e-02 | +1.728e-01 | 1.039e-04 | 7.652e-10 | -2.35e-09 |  |
| pair | (1,+1) -> (5,-21) + (6,-9) | -4.449e-01 | +6.430e-01 | 1.250e-04 | 4.939e-09 | -2.34e-09 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | runaway | 7.09e+09 | 6.0863e-04 | 15.4 | 0.002531 | 2 |
| d3 | 3 | 2 |  |  | runaway | 8.09e+09 | 5.4522e-04 | 16.2 | 0.002698 | 3 |
| d5 | 5 | 3 |  |  | runaway | 7.93e+09 | 5.2497e-04 | 17.9 | 0.001321 | 5 |
| d6 | 6 | 4 |  |  | runaway | 8.74e+09 | 4.7993e-04 | 16.5 | 0.00461 | 6 |
| d7 | 7 | 5 |  |  | runaway | 8.76e+09 | 5.0060e-04 | 18.6 | 0.002749 | 5 |
| d8 | 8 | 6 |  |  | runaway | 9.41e+09 | 2.5824e-04 | 18.2 | 0.003507 | 8 |
| d10 | 10 | 7 |  |  | runaway | 1.04e+10 | 4.1768e-04 | 10.7 | 0.01538 | 10 |
| d11 | 11 | 8 |  |  | runaway | 1.04e+10 | 2.7434e-04 | 8.48 | 0.01006 | 11 |
| d12 | 12 | 9 |  |  | runaway | 1.01e+10 | 4.2980e-04 | 12.5 | 0.01018 | 12 |
| d13 | 13 | 10 |  |  | runaway | 1.19e+10 | 3.5222e-04 | 16.1 | 0.01086 | 13 |
| d15 | 15 | 11 |  |  | runaway | 1.17e+10 | 3.6405e-04 | 14.5 | 0.01049 | 15 |
| d17 | 17 | 12 |  |  | runaway | 1.34e+10 | 2.9874e-04 | 16 | 0.01086 | 17 |
| d18 | 18 | 13 |  |  | runaway | 1.3e+10 | 3.3080e-04 | 12.2 | 0.01064 | 18 |
| d19 | 19 | 14 |  |  | runaway | 1.27e+10 | 3.0415e-05 | 9.9 | 0.01811 | 19 |
| d21 | 21 | 15 |  |  | runaway | 1.39e+10 | 3.0851e-04 | 16.9 | 0.01007 | 21 |
| d23 | 23 | 16 |  |  | runaway | 1.47e+10 | 2.9406e-04 | 13.5 | 0.0105 | 23 |
| d24 | 24 | 17 |  |  | runaway | 1.43e+10 | 3.0453e-04 | 14.8 | 0.009872 | 24 |

Same run from some N on: not reached.

### ladder_1p/l1_n+2

Parents: (1,+2) f 16.535 c/d, gamma -2.622e-09 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (1,+2) -> (1,-2) + (2,-4) | +2.919e-03 | -1.322e-01 | 3.898e-07 | 1.386e-14 | -2.62e-09 |  |
| pair | (1,+2) -> (3,-6) + (4,-7) | -5.871e-03 | +3.878e-02 | 2.080e-05 | 4.717e-12 | -2.62e-09 |  |
| pair | (1,+2) -> (5,-9) + (6,-11) | +4.997e-02 | -1.595e-01 | 8.024e-05 | 1.240e-10 | -2.62e-09 |  |
| pair | (1,+2) -> (5,-10) + (6,-10) | -1.707e-03 | +4.775e-03 | 1.174e-04 | 1.745e-10 | -2.62e-09 |  |
| pair | (1,+2) -> (5,-8) + (6,-12) | +1.177e-01 | -2.504e-01 | 1.480e-04 | 1.987e-10 | -2.62e-09 |  |
| pair | (1,+2) -> (4,-9) + (5,-8) | -1.819e-01 | -2.502e-01 | 4.835e-04 | 2.327e-10 | -2.62e-09 |  |
| pair | (1,+2) -> (5,-8) + (6,-13) | -2.031e-01 | +2.263e-01 | 5.172e-04 | 8.866e-10 | -2.62e-09 |  |
| pair | (1,+2) -> (5,-11) + (6,-9) | +2.220e-01 | +2.012e-01 | 7.928e-04 | 8.260e-10 | -2.62e-09 |  |
| pair | (1,+2) -> (5,-12) + (6,-9) | -3.234e-01 | -1.882e-01 | 1.920e-03 | 2.520e-09 | -2.62e-09 |  |
| pair | (1,+2) -> (2,-4) + (3,-6) | -3.453e-01 | -1.551e-01 | 4.417e-03 | 5.202e-10 | -2.62e-09 |  |
| pair | (1,+2) -> (4,-8) + (5,-9) | -3.233e-01 | +6.671e-02 | 2.215e-02 | 1.296e-08 | -2.62e-09 |  |
| pair | (1,+2) -> (4,-8) + (5,-8) | +3.042e-01 | +4.602e-02 | 3.848e-02 | 1.491e-08 | -2.62e-09 |  |
| pair | (1,+2) -> (4,-7) + (5,-10) | -3.053e-01 | -2.567e-02 | 1.136e-01 | 8.442e-08 | -2.62e-09 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | runaway | 5.81e+09 | 6.6463e-04 | 14.5 | 4.737e-05 | 2 |
| d4 | 4 | 2 |  |  | runaway | 6.53e+09 | 2.2490e-04 | 11.8 | 0.0003855 | 4 |
| d6 | 6 | 3 |  |  | runaway | 7.21e+09 | 5.3318e-04 | 13.6 | 0.0005746 | 6 |
| d8 | 8 | 4 |  |  | runaway | 7.1e+09 | 4.7635e-04 | 11.1 | 0.0006134 | 6 |
| d10 | 10 | 5 |  |  | runaway | 6.74e+09 | 5.6532e-04 | 8.91 | 0.0007058 | 8 |
| d11 | 11 | 6 |  |  | runaway | 8.45e+09 | 4.5584e-04 | 13.2 | 0.001452 | 9 |
| d12 | 12 | 7 |  |  | runaway | 8.46e+09 | 4.4172e-04 | 9.39 | 0.001787 | 10 |
| d14 | 14 | 8 |  |  | runaway | 8.58e+09 | 4.3382e-04 | 10.9 | 0.001748 | 12 |
| d15 | 15 | 10 |  |  | runaway | 9.62e+09 | 3.9940e-04 | 12.1 | 0.002522 | 13 |
| d16 | 16 | 13 |  |  | runaway | 1.13e+10 | 3.0638e-04 | 9.86 | 0.004834 | 16 |

Same run from some N on: not reached.

### ladder_1p/l1_n+3

Parents: (1,+3) f 18.061 c/d, gamma -1.983e-08 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (1,+3) -> (4,-7) + (5,-8) | -3.660e-01 | -1.300e-01 | 6.331e-03 | 3.247e-11 | -1.98e-08 |  |
| pair | (1,+3) -> (5,-9) + (6,-9) | +2.518e-01 | -5.939e-02 | 1.237e-02 | 9.409e-11 | -1.98e-08 |  |
| pair | (1,+3) -> (5,-8) + (6,-10) | -6.247e-02 | -6.130e-03 | 6.896e-02 | 9.724e-10 | -1.98e-08 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | runaway | 5.82e+08 | 8.7341e-04 | 10.1 | 1.538e-06 | 2 |
| d4 | 4 | 2 |  |  | runaway | 5.82e+08 | 8.7341e-04 | 10.1 | 1.538e-06 | 2 |
| d5 | 5 | 3 |  |  | runaway | 5.82e+08 | 8.7342e-04 | 10.1 | 1.589e-06 | 3 |

Same run from some N on: not reached.

### ladder_1p/l1_n+7

Parents: (1,+7) f 32.096 c/d, gamma -7.808e-07 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| harmonic | (1,+7) + (1,+7) -> (0,+16) | +7.401e-03 | +1.147e+02 |  | 6.687e-11 | +2.30e-05 | 2.181e+04 |
| harmonic | (1,+7) + (1,+7) -> (2,+15) | -2.827e-01 | +1.103e+02 |  | 1.305e-10 | +2.27e-05 | 1.612e+04 |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  | (1,+7) + (1,+7) -> (0,+16) | end | 1e+13 | 6.6875e-11 | 3.65e-10 | 1 | 1 |
| d2 | 2 | 2 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (1,+7) -> (2,+15) | end | 1e+13 | 4.3618e-11 | 8.64e-10 | 0.9985 | 2 |

Same run from some N on: not reached.

### ladder_1p/l1_n+8

Parents: (1,+8) f 35.663 c/d, gamma -3.347e-10 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| harmonic | (1,+8) + (1,+8) -> (2,+17) | +1.731e-01 | +9.044e+01 |  | 5.396e-14 | +3.03e-05 | 1.432e+04 |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  | (1,+8) + (1,+8) -> (2,+17) | end | 1e+13 | 5.3960e-14 | 9.33e-12 | 1.002 | 1 |

Same run from some N on: not reached.

### ladder_1p/l1_n-1

Parents: (1,-1) f 10.979 c/d, gamma -1.166e-10 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (1,-1) -> (5,-8) + (6,-58) | +3.229e-03 | +5.637e-02 | 5.113e-09 | 5.285e-09 | +7.02e-09 |  |
| pair | (1,-1) -> (5,-10) + (6,-31) | +7.958e-04 | +1.392e-01 | 7.018e-09 | 1.719e-08 | -4.55e-11 |  |
| pair | (1,-1) -> (4,-29) + (5,-9) | -7.595e-03 | +3.436e-01 | 7.368e-08 | 6.770e-07 | -2.89e-11 |  |
| pair | (1,-1) -> (4,-40) + (5,-8) | -1.035e-02 | +1.890e-01 | 7.511e-08 | 1.375e-07 | +3.30e-10 |  |
| pair | (1,-1) -> (6,-10) + (7,-56) | -5.880e-03 | +4.826e-02 | 8.064e-08 | 8.565e-08 | +3.81e-09 |  |
| pair | (1,-1) -> (6,-14) + (7,-27) | +2.035e-03 | +1.458e-01 | 1.322e-07 | 5.713e-08 | -7.04e-11 |  |
| pair | (1,-1) -> (4,-7) + (5,-35) | +1.108e-02 | +1.884e-01 | 1.622e-07 | 1.554e-06 | +5.57e-11 |  |
| pair | (1,-1) -> (5,-41) + (6,-10) | -1.846e-02 | +1.979e-01 | 2.930e-07 | 4.398e-07 | +5.18e-10 |  |
| pair | (1,-1) -> (5,-29) + (6,-11) | +9.849e-03 | +2.761e-01 | 3.305e-07 | 2.363e-06 | -3.17e-11 |  |
| pair | (1,-1) -> (4,-39) + (5,-8) | +2.224e-02 | +1.942e-01 | 3.460e-07 | 6.684e-07 | +2.99e-10 |  |
| pair | (1,-1) -> (5,-9) + (6,-42) | +2.257e-02 | +1.326e-01 | 5.528e-07 | 7.988e-07 | +5.76e-10 |  |
| pair | (1,-1) -> (5,-8) + (6,-57) | +3.643e-02 | +5.469e-02 | 9.134e-07 | 9.549e-07 | +5.18e-09 |  |
| pair | (1,-1) -> (6,-10) + (7,-55) | +2.893e-02 | +5.941e-02 | 1.031e-06 | 1.082e-06 | +4.71e-09 |  |
| pair | (1,-1) -> (4,-41) + (5,-8) | -5.270e-02 | +2.118e-01 | 1.039e-06 | 1.512e-06 | +5.66e-10 |  |
| pair | (1,-1) -> (4,-7) + (5,-36) | -3.857e-02 | +2.055e-01 | 1.180e-06 | 4.253e-06 | +1.30e-10 |  |
| pair | (1,-1) -> (3,-22) + (4,-7) | -1.838e-02 | +4.149e-01 | 1.189e-06 | 1.819e-08 | -1.04e-10 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | end | 1e+13 | 6.2713e-09 | 752 | 9.601e+281 | 1 |
| d4 | 4 | 2 |  |  | end | 1e+13 | 6.3828e-09 | 72.9 | 0.8958 | 3 |
| d6 | 6 | 3 |  |  | end | 1e+13 | 5.5360e-09 | 5.38 | 0.906 | 2 |
| d7 | 7 | 4 |  |  | end | 1e+13 | 4.5062e-09 | 4.51 | 1.123 | 2 |
| d9 | 9 | 5 |  |  | end | 1e+13 | 4.5062e-09 | 4.51 | 1.123 | 2 |
| d11 | 11 | 6 |  |  | end | 1e+13 | 4.7034e-09 | 4.15 | 0.9453 | 2 |
| d13 | 13 | 7 |  |  | end | 1e+13 | 4.7034e-09 | 4.15 | 0.9453 | 2 |
| d14 | 14 | 8 |  |  | end | 1e+13 | 4.7033e-09 | 4.15 | 0.9453 | 2 |
| d16 | 16 | 9 |  |  | end | 1e+13 | 4.7034e-09 | 4.15 | 0.9453 | 2 |
| d17 | 17 | 10 |  |  | end | 1e+13 | 4.6444e-09 | 5.53 | 1.082 | 2 |
| d18 | 18 | 11 |  |  | end | 1e+13 | 4.6851e-09 | 5.47 | 1.035 | 2 |
| d19 | 19 | 12 |  |  | end | 1e+13 | 4.5863e-09 | 5.69 | 1.165 | 2 |
| d20 | 20 | 13 |  |  | end | 1e+13 | 4.5863e-09 | 5.69 | 1.165 | 2 |
| d21 | 21 | 14 |  |  | end | 1e+13 | 4.6681e-09 | 5.01 | 0.8967 | 2 |
| d22 | 22 | 15 |  |  | end | 1e+13 | 4.6681e-09 | 5.01 | 0.8967 | 2 |
| d23 | 23 | 16 |  |  | end | 1e+13 | 4.6681e-09 | 5.01 | 0.8967 | 2 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 21 on.

### ladder_1p/l2_n+0

Parents: (2,+0) f 14.909 c/d, gamma -2.632e-09 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (2,+0) -> (6,-19) + (6,-9) | -1.729e-02 | +2.204e-01 | 1.600e-06 | 2.943e-11 | -2.62e-09 |  |
| pair | (2,+0) -> (6,-9) + (8,-26) | -2.043e-01 | +4.616e-01 | 1.991e-05 | 2.816e-09 | -2.60e-09 |  |
| pair | (2,+0) -> (6,-9) + (8,-27) | -3.336e-01 | +5.949e-01 | 2.124e-05 | 7.376e-09 | -2.58e-09 |  |
| pair | (2,+0) -> (6,-12) + (6,-12) | +4.555e-02 | +1.969e-01 | 5.984e-05 | 1.831e-10 | -2.63e-09 |  |
| pair | (2,+0) -> (3,-7) + (3,-6) | -4.770e-02 | -1.202e-01 | 1.766e-04 | 4.283e-11 | -2.63e-09 |  |
| pair | (2,+0) -> (4,-7) + (6,-15) | -1.371e-01 | +2.026e-01 | 2.205e-04 | 1.109e-09 | -2.63e-09 |  |
| pair | (2,+0) -> (6,-20) + (6,-9) | -2.873e-01 | +2.596e-01 | 2.756e-04 | 7.718e-09 | -2.62e-09 |  |
| pair | (2,+0) -> (4,-10) + (4,-7) | +9.317e-02 | +1.093e-01 | 7.122e-04 | 4.780e-10 | -2.63e-09 |  |
| pair | (2,+0) -> (5,-12) + (5,-9) | +8.565e-02 | +9.938e-02 | 7.307e-04 | 1.289e-09 | -2.63e-09 |  |
| pair | (2,+0) -> (6,-18) + (6,-9) | +2.593e-01 | +1.727e-01 | 7.438e-04 | 7.304e-09 | -2.62e-09 |  |
| pair | (2,+0) -> (4,-8) + (6,-12) | +2.321e-01 | +2.290e-01 | 9.038e-04 | 1.334e-09 | -2.63e-09 |  |
| pair | (2,+0) -> (6,-14) + (6,-11) | -1.697e-01 | +1.830e-01 | 9.112e-04 | 4.891e-09 | -2.63e-09 |  |
| pair | (2,+0) -> (4,-9) + (4,-8) | -6.747e-02 | -7.352e-02 | 9.484e-04 | 5.352e-10 | -2.63e-09 |  |
| pair | (2,+0) -> (3,-6) + (5,-11) | +9.265e-02 | +7.259e-02 | 1.230e-03 | 1.221e-09 | -2.63e-09 |  |
| pair | (2,+0) -> (5,-14) + (5,-8) | +5.791e-02 | +4.229e-02 | 1.238e-03 | 3.414e-09 | -2.63e-09 |  |
| pair | (2,+0) -> (4,-10) + (6,-11) | -3.893e-01 | +3.754e-01 | 1.245e-03 | 2.358e-09 | -2.63e-09 |  |
| pair | (2,+0) -> (4,-12) + (6,-9) | +2.239e-01 | +1.741e-01 | 1.311e-03 | 1.642e-09 | -2.63e-09 |  |
| pair | (2,+0) -> (4,-13) + (6,-9) | -5.945e-02 | -4.716e-02 | 1.341e-03 | 1.640e-09 | -2.63e-09 |  |
| pair | (2,+0) -> (6,-15) + (6,-10) | +1.665e-01 | +1.191e-01 | 1.600e-03 | 1.082e-08 | -2.62e-09 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | runaway | 5.87e+09 | 6.6573e-04 | 11.5 | 0.001852 | 2 |
| d3 | 3 | 2 |  |  | runaway | 5.97e+09 | 6.4316e-04 | 5.85 | 0.004709 | 3 |
| d4 | 4 | 3 |  |  | runaway | 5.94e+09 | 6.5130e-04 | 5.8 | 0.006605 | 4 |
| d5 | 5 | 4 |  |  | runaway | 6.87e+09 | 5.6448e-04 | 8.53 | 0.006429 | 5 |
| d7 | 7 | 5 |  |  | runaway | 6.79e+09 | 5.7235e-04 | 9 | 0.00618 | 5 |
| d9 | 9 | 6 |  |  | runaway | 7.01e+09 | 5.4386e-04 | 8.71 | 0.006553 | 7 |
| d10 | 10 | 7 |  |  | runaway | 6.27e+09 | 6.1742e-04 | 8.05 | 0.006121 | 6 |
| d11 | 11 | 8 |  |  | runaway | 6.71e+09 | 5.7174e-04 | 8.51 | 0.006031 | 11 |
| d13 | 13 | 9 |  |  | runaway | 7.08e+09 | 5.4661e-04 | 7.46 | 0.006342 | 13 |
| d14 | 14 | 10 |  |  | runaway | 8.64e+09 | 3.2762e-04 | 11.9 | 0.006396 | 14 |
| d15 | 15 | 11 |  |  | runaway | 1.1e+10 | 3.3802e-04 | 10.4 | 0.006963 | 15 |
| d17 | 17 | 12 |  |  | runaway | 1.14e+10 | 3.2421e-04 | 10.5 | 0.007196 | 17 |
| d18 | 18 | 13 |  |  | runaway | 1.13e+10 | 2.7314e-04 | 10.2 | 0.007858 | 18 |
| d19 | 19 | 14 |  |  | runaway | 1.19e+10 | 3.1002e-04 | 15.6 | 0.006646 | 19 |
| d21 | 21 | 16 |  |  | runaway | 1.35e+10 | 2.2024e-04 | 12.8 | 0.008181 | 21 |
| d22 | 22 | 17 |  |  | runaway | 1.3e+10 | 2.8194e-04 | 10.5 | 0.006727 | 22 |
| d23 | 23 | 18 |  |  | runaway | 1.28e+10 | 2.8313e-04 | 13.2 | 0.005993 | 23 |
| d24 | 24 | 19 |  |  | runaway | 1.83e+10 | 1.9314e-04 | 13.1 | 0.01577 | 24 |

Same run from some N on: not reached.

### ladder_1p/l2_n+1

Parents: (2,+1) f 17.181 c/d, gamma -1.006e-08 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (2,+1) -> (3,-6) + (5,-8) | +3.162e-03 | -1.312e-01 | 4.852e-07 | 8.150e-15 | -1.01e-08 |  |
| pair | (2,+1) -> (5,-9) + (5,-9) | -8.600e-02 | -1.054e-01 | 5.699e-04 | 2.775e-11 | -1.01e-08 |  |
| pair | (2,+1) -> (4,-7) + (4,-7) | -1.427e-01 | -1.155e-01 | 1.314e-03 | 2.428e-11 | -1.01e-08 |  |
| pair | (2,+1) -> (4,-9) + (6,-9) | -2.392e-01 | +6.755e-02 | 1.036e-02 | 2.960e-10 | -1.01e-08 |  |
| pair | (2,+1) -> (6,-12) + (6,-9) | +6.036e-02 | -1.306e-02 | 1.223e-02 | 1.030e-09 | -1.01e-08 |  |
| pair | (2,+1) -> (6,-13) + (6,-9) | -2.604e-01 | +5.220e-02 | 1.353e-02 | 1.468e-09 | -1.01e-08 |  |
| pair | (2,+1) -> (4,-8) + (6,-9) | +2.469e-01 | -5.878e-02 | 1.425e-02 | 3.229e-10 | -1.01e-08 |  |
| pair | (2,+1) -> (6,-11) + (6,-10) | -3.216e-01 | +7.605e-02 | 1.545e-02 | 2.165e-09 | -1.01e-08 |  |
| pair | (2,+1) -> (4,-7) + (6,-10) | +1.609e-01 | -2.318e-02 | 3.346e-02 | 1.752e-09 | -1.01e-08 |  |
| pair | (2,+1) -> (5,-10) + (5,-8) | -2.962e-01 | +2.505e-02 | 1.073e-01 | 5.675e-09 | -1.01e-08 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | runaway | 1.41e+09 | 7.0500e-04 | 14.1 | 2.988e-05 | 2 |
| d3 | 3 | 2 |  |  | runaway | 1.54e+09 | 6.3645e-04 | 8.61 | 0.000163 | 3 |
| d4 | 4 | 3 |  |  | runaway | 1.54e+09 | 5.9116e-04 | 8.27 | 0.0001524 | 4 |
| d6 | 6 | 4 |  |  | runaway | 1.54e+09 | 6.4601e-04 | 9.42 | 0.0001202 | 6 |
| d7 | 7 | 5 |  |  | runaway | 1.56e+09 | 6.4722e-04 | 7.04 | 0.0001691 | 5 |
| d8 | 8 | 6 |  |  | runaway | 1.53e+09 | 6.5383e-04 | 8.7 | 0.0001315 | 8 |
| d9 | 9 | 7 |  |  | runaway | 1.59e+09 | 6.3188e-04 | 9.22 | 0.0002328 | 9 |
| d11 | 11 | 9 |  |  | runaway | 1.74e+09 | 4.2226e-04 | 9.9 | 0.000545 | 11 |
| d12 | 12 | 10 |  |  | runaway | 1.74e+09 | 5.7690e-04 | 8.76 | 0.0007628 | 12 |

Same run from some N on: not reached.

### ladder_1p/l2_n+2

Parents: (2,+2) f 19.961 c/d, gamma -4.601e-08 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| harmonic | (2,+2) + (2,+2) -> (4,+8) | -1.189e-01 | +9.158e+00 |  | 2.607e-09 | +2.07e-06 | 2.975e+03 |
| pair | (2,+2) -> (6,-9) + (6,-9) | -4.321e-01 | +5.133e-02 | 4.644e-02 | 3.363e-11 | -4.60e-08 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  | (2,+2) + (2,+2) -> (4,+8) | end | 1e+13 | 2.6074e-09 | 1.13e-10 | 0.997 | 1 |
| d2 | 2 | 2 |  | (2,+2) + (2,+2) -> (4,+8) | end | 1e+13 | 2.6074e-09 | 1.11e-10 | 0.997 | 1 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 1 on.

### ladder_1p/l2_n+5

Parents: (2,+5) f 28.190 c/d, gamma -8.234e-08 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| harmonic | (2,+5) + (2,+5) -> (0,+14) | +1.888e-01 | +1.118e+01 |  | 1.108e-09 | +1.70e-05 | 2.089e+03 |
| harmonic | (2,+5) + (2,+5) -> (2,+13) | -1.047e-01 | +7.492e+00 |  | 1.785e-09 | +1.67e-05 | 1.653e+03 |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  | (2,+5) + (2,+5) -> (0,+14) | end | 1e+13 | 1.1083e-09 | 5.44e-11 | 1.003 | 1 |
| d2 | 2 | 2 |  | (2,+5) + (2,+5) -> (0,+14); (2,+5) + (2,+5) -> (2,+13) | end | 1e+13 | 6.8030e-10 | 3.42e-11 | 1.001 | 2 |

Same run from some N on: not reached.

### ladder_1p/l2_n+6

Parents: (2,+6) f 30.298 c/d, gamma -7.891e-07 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| harmonic | (2,+6) + (2,+6) -> (0,+15) | -2.039e-01 | +1.112e+02 |  | 1.096e-10 | +1.95e-05 | 1.897e+04 |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  | (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 1.0963e-10 | 6.75e-11 | 0.9966 | 1 |

Same run from some N on: not reached.

### ladder_1p/l2_n+7

Parents: (2,+7) f 33.770 c/d, gamma -5.357e-07 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| harmonic | (2,+7) + (2,+7) -> (2,+16) | +1.674e-01 | +6.838e+01 |  | 1.577e-10 | +2.64e-05 | 1.122e+04 |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  | (2,+7) + (2,+7) -> (2,+16) | end | 1e+13 | 1.5771e-10 | 9.39e-10 | 1.002 | 1 |

Same run from some N on: not reached.

### ladder_1p/l2_n-1

Parents: (2,-1) f 12.733 c/d, gamma -4.582e-10 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (2,-1) -> (6,-9) + (8,-45) | +3.115e-02 | +2.976e-01 | 5.920e-08 | 2.098e-07 | +5.19e-10 |  |
| pair | (2,-1) -> (6,-34) + (6,-9) | +1.688e-02 | +3.345e-01 | 6.888e-08 | 3.827e-08 | -2.63e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-46) | -3.960e-02 | +2.907e-01 | 8.054e-08 | 2.013e-07 | +7.88e-10 |  |
| pair | (2,-1) -> (4,-10) + (6,-15) | +1.920e-03 | +1.635e-01 | 1.461e-07 | 3.195e-11 | -4.52e-10 |  |
| pair | (2,-1) -> (5,-8) + (7,-33) | -1.907e-02 | +2.829e-01 | 1.607e-07 | 4.596e-08 | -2.99e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-47) | -1.035e-01 | +2.555e-01 | 7.101e-07 | 1.726e-06 | +8.20e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-44) | +9.324e-02 | +2.573e-01 | 1.162e-06 | 2.499e-05 | +1.26e-10 |  |
| pair | (2,-1) -> (4,-24) + (6,-9) | -6.159e-02 | +7.050e-01 | 1.419e-06 | 6.334e-09 | -4.30e-10 |  |
| pair | (2,-1) -> (6,-35) + (6,-9) | -5.708e-02 | +2.698e-01 | 1.555e-06 | 4.151e-07 | -3.02e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-49) | -1.989e-01 | +2.514e-01 | 1.843e-06 | 3.156e-06 | +1.48e-09 |  |
| pair | (2,-1) -> (5,-24) + (5,-8) | -2.537e-02 | +2.504e-01 | 1.944e-06 | 9.044e-09 | -4.29e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-43) | +1.437e-01 | +2.842e-01 | 2.174e-06 | 3.953e-05 | +1.40e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-48) | -1.489e-01 | +2.143e-01 | 2.366e-06 | 6.557e-06 | +6.89e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-50) | -2.587e-01 | +2.448e-01 | 2.584e-06 | 3.854e-06 | +2.07e-09 |  |
| pair | (2,-1) -> (6,-10) + (8,-34) | -4.684e-02 | +2.261e-01 | 3.028e-06 | 9.865e-07 | -2.92e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-42) | +2.166e-01 | +3.337e-01 | 3.187e-06 | 3.463e-05 | +1.99e-10 |  |
| pair | (2,-1) -> (6,-33) + (6,-9) | +1.151e-01 | +3.510e-01 | 3.249e-06 | 1.118e-06 | -2.89e-10 |  |
| pair | (2,-1) -> (4,-23) + (6,-9) | +6.999e-02 | +6.117e-01 | 3.412e-06 | 6.649e-09 | -4.39e-10 |  |
| pair | (2,-1) -> (6,-36) + (6,-9) | -1.161e-01 | +3.067e-01 | 3.424e-06 | 3.596e-06 | -2.26e-10 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d4 | 4 | 3 |  |  | end | 1e+13 | 2.6766e-08 | 12.2 | 0.9845 | 1 |
| d6 | 6 | 4 |  |  | end | 1e+13 | 2.7354e-08 | 5.12 | 1.004 | 1 |
| d8 | 8 | 5 |  |  | end | 1e+13 | 2.6845e-08 | 2.69 | 0.9948 | 2 |
| d9 | 9 | 6 |  |  | end | 1e+13 | 2.6143e-08 | 3.03 | 0.9983 | 1 |
| d10 | 10 | 7 |  |  | end | 1e+13 | 2.5762e-08 | 3.03 | 0.9992 | 2 |
| d11 | 11 | 8 |  |  | end | 1e+13 | 2.0524e-08 | 9.18 | 0.9903 | 1 |
| d12 | 12 | 9 |  |  | end | 1e+13 | 1.9758e-08 | 5.39 | 0.993 | 2 |
| d13 | 13 | 10 |  |  | end | 1e+13 | 2.0069e-08 | 2.98 | 0.9983 | 2 |
| d14 | 14 | 11 |  |  | end | 1e+13 | 1.9680e-08 | 0.357 | 1 | 2 |
| d15 | 15 | 12 |  |  | end | 1e+13 | 1.9761e-08 | 0.427 | 0.9997 | 1 |
| d16 | 16 | 13 |  |  | end | 1e+13 | 1.9282e-08 | 2.18 | 0.9973 | 1 |
| d17 | 17 | 14 |  |  | end | 1e+13 | 1.9341e-08 | 0.988 | 0.9999 | 2 |
| d19 | 19 | 15 |  |  | end | 1e+13 | 1.9130e-08 | 0.798 | 1 | 2 |
| d20 | 20 | 16 |  |  | end | 1e+13 | 1.9260e-08 | 1.34 | 0.9996 | 1 |
| d21 | 21 | 17 |  |  | end | 1e+13 | 2.0269e-08 | 4.15 | 0.9965 | 1 |
| d22 | 22 | 18 |  |  | end | 1e+13 | 2.2537e-08 | 0.358 | 0.9997 | 1 |
| d23 | 23 | 19 |  |  | end | 1e+13 | 2.2117e-08 | 0.468 | 1 | 1 |

Same run from some N on: not reached.

### ladder_1p/l2_n-1_ladder2

Parents: (2,-1) f 12.733 c/d, gamma -4.582e-10 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first; only pairs with an E_eq.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (2,-1) -> (6,-9) + (8,-45) | +3.115e-02 | +2.976e-01 | 5.920e-08 | 2.098e-07 | +5.19e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-46) | -3.960e-02 | +2.907e-01 | 8.054e-08 | 2.013e-07 | +7.88e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-47) | -1.035e-01 | +2.555e-01 | 7.101e-07 | 1.726e-06 | +8.20e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-44) | +9.324e-02 | +2.573e-01 | 1.162e-06 | 2.499e-05 | +1.26e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-49) | -1.989e-01 | +2.514e-01 | 1.843e-06 | 3.156e-06 | +1.48e-09 |  |
| pair | (2,-1) -> (6,-9) + (8,-43) | +1.437e-01 | +2.842e-01 | 2.174e-06 | 3.953e-05 | +1.40e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-48) | -1.489e-01 | +2.143e-01 | 2.366e-06 | 6.557e-06 | +6.89e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-50) | -2.587e-01 | +2.448e-01 | 2.584e-06 | 3.854e-06 | +2.07e-09 |  |
| pair | (2,-1) -> (6,-9) + (8,-42) | +2.166e-01 | +3.337e-01 | 3.187e-06 | 3.463e-05 | +1.99e-10 |  |
| pair | (2,-1) -> (6,-9) + (8,-51) | -3.149e-01 | +2.207e-01 | 4.573e-06 | 6.676e-06 | +2.20e-09 |  |
| pair | (2,-1) -> (6,-9) + (8,-41) | +3.009e-01 | +3.461e-01 | 6.580e-06 | 2.120e-04 | +9.80e-11 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d11 | 11 | 10 |  |  | end | 1e+13 | 3.0790e-08 | 51.6 | 0.889 | 1 |
| d12 | 12 | 11 |  |  | end | 1e+13 | 2.9188e-08 | 35.2 | 1.146 | 1 |

Same run from some N on: not reached.

### ladder_1p/l2_n-2

Parents: (2,-2) f 10.944 c/d, gamma -6.697e-11 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (2,-2) -> (3,-31) + (5,-8) | -4.038e-05 | +2.147e-01 | 4.870e-12 | 1.499e-10 | +1.47e-11 |  |
| pair | (2,-2) -> (5,-9) + (7,-50) | +8.218e-04 | +5.254e-02 | 1.302e-09 | 1.373e-09 | +2.48e-09 |  |
| pair | (2,-2) -> (4,-29) + (4,-7) | -1.662e-03 | +1.754e-01 | 8.448e-09 | 1.563e-07 | +2.03e-11 |  |
| pair | (2,-2) -> (4,-24) + (6,-11) | -1.126e-03 | +2.868e-01 | 1.088e-08 | 7.362e-09 | -3.67e-11 |  |
| pair | (2,-2) -> (5,-21) + (5,-11) | -8.884e-04 | +1.051e-01 | 7.174e-08 | 7.245e-09 | -5.08e-11 |  |
| pair | (2,-2) -> (6,-29) + (6,-12) | +2.325e-03 | +8.988e-02 | 1.762e-07 | 5.422e-06 | +1.47e-11 |  |
| pair | (2,-2) -> (6,-49) + (6,-10) | +7.677e-03 | +5.037e-02 | 2.559e-07 | 2.744e-07 | +1.89e-09 |  |
| pair | (2,-2) -> (6,-35) + (6,-11) | +3.389e-03 | +6.764e-02 | 3.602e-07 | 1.089e-06 | +9.07e-11 |  |
| pair | (2,-2) -> (3,-6) + (5,-29) | +1.660e-02 | +1.748e-01 | 6.420e-07 | 1.674e-05 | +1.63e-11 |  |
| pair | (2,-2) -> (4,-13) + (4,-11) | -2.590e-03 | +1.458e-01 | 6.569e-07 | 2.949e-09 | -6.28e-11 |  |
| pair | (2,-2) -> (6,-10) + (8,-64) | +1.326e-02 | +2.058e-02 | 7.872e-07 | 7.966e-07 | +1.13e-08 |  |
| pair | (2,-2) -> (2,-4) + (4,-25) | +1.402e-02 | +1.489e-01 | 8.751e-07 | 1.326e-06 | -3.00e-11 |  |
| pair | (2,-2) -> (4,-7) + (6,-42) | +2.850e-02 | +1.035e-01 | 8.953e-07 | 1.096e-06 | +6.25e-10 |  |
| pair | (2,-2) -> (6,-10) + (8,-65) | -1.328e-02 | +1.947e-02 | 9.589e-07 | 9.708e-07 | +1.05e-08 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | end | 1e+13 | 1.4977e-10 | 0.00843 | 1 | 2 |
| d4 | 4 | 2 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d6 | 6 | 3 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d8 | 8 | 4 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d10 | 10 | 5 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d12 | 12 | 6 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d14 | 14 | 7 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d15 | 15 | 8 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d16 | 16 | 9 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d17 | 17 | 9 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d19 | 19 | 10 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d20 | 20 | 11 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d22 | 22 | 12 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d23 | 23 | 13 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d24 | 24 | 14 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 2 on.

### ladder_1p/l2_n-3

Parents: (2,-3) f 9.373 c/d, gamma -8.642e-13 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (2,-3) -> (14,-48) + (14,-48) | -5.624e-05 | -4.521e-01 | 6.397e-11 | 6.400e-11 | +2.75e-09 |  |
| pair | (2,-3) -> (14,-44) + (16,-60) | -4.809e-06 | -9.464e-02 | 1.092e-09 | 1.092e-09 | +6.30e-09 |  |
| pair | (2,-3) -> (11,-39) + (13,-43) | -4.959e-04 | -6.603e-01 | 1.419e-09 | 1.423e-09 | +7.29e-10 |  |
| pair | (2,-3) -> (2,-17) + (4,-9) | +3.922e-04 | +3.763e-01 | 2.489e-09 | 3.364e-09 | +5.31e-12 |  |
| pair | (2,-3) -> (9,-32) + (9,-30) | -3.111e-04 | -2.937e-01 | 2.923e-09 | 2.950e-09 | +1.85e-10 |  |
| pair | (2,-3) -> (3,-13) + (5,-15) | +5.393e-04 | +2.764e-01 | 9.578e-09 | 1.249e-08 | +6.08e-12 |  |
| pair | (2,-3) -> (16,-54) + (18,-62) | +5.794e-04 | -2.621e-01 | 1.261e-08 | 1.261e-08 | +1.37e-08 |  |
| pair | (2,-3) -> (12,-44) + (14,-45) | -6.520e-04 | -2.808e-01 | 1.476e-08 | 1.477e-08 | +1.80e-09 |  |
| pair | (2,-3) -> (6,-17) + (8,-35) | +9.380e-04 | -1.569e-01 | 1.658e-08 | 1.680e-08 | +1.35e-10 |  |
| pair | (2,-3) -> (9,-33) + (9,-29) | -5.967e-04 | +1.702e-01 | 3.179e-08 | 3.205e-08 | +2.26e-10 |  |
| pair | (2,-3) -> (16,-55) + (18,-61) | -1.486e-03 | -3.945e-01 | 3.438e-08 | 3.440e-08 | +1.18e-08 |  |
| pair | (2,-3) -> (4,-26) + (6,-14) | +7.398e-04 | +1.153e-01 | 4.925e-08 | 5.140e-08 | +4.07e-11 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  |  | end | 1e+13 | 6.4083e-11 | 0.196 | 1.009 | 1 |
| d3 | 3 | 2 |  |  | end | 1e+13 | 6.4014e-11 | 0.0729 | 0.9984 | 1 |
| d5 | 5 | 3 |  |  | end | 1e+13 | 6.4057e-11 | 0.173 | 1.006 | 1 |
| d7 | 7 | 4 |  |  | end | 1e+13 | 6.4057e-11 | 0.173 | 1.006 | 1 |
| d9 | 9 | 5 |  |  | end | 1e+13 | 6.4057e-11 | 0.173 | 1.006 | 1 |
| d11 | 11 | 6 |  |  | end | 1e+13 | 6.4057e-11 | 0.173 | 1.006 | 1 |
| d13 | 13 | 7 |  |  | end | 1e+13 | 6.4057e-11 | 0.173 | 1.006 | 1 |
| d15 | 15 | 8 |  |  | end | 1e+13 | 6.4057e-11 | 0.173 | 1.006 | 1 |
| d17 | 17 | 9 |  |  | end | 1e+13 | 6.4057e-11 | 0.173 | 1.006 | 1 |
| d19 | 19 | 10 |  |  | end | 1e+13 | 6.4057e-11 | 0.173 | 1.006 | 1 |
| d21 | 21 | 13 |  |  | end | 1e+13 | 6.4057e-11 | 0.173 | 1.006 | 1 |
| d23 | 23 | 14 |  |  | end | 1e+13 | 6.4057e-11 | 0.173 | 1.006 | 1 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 5 on.

## ladder_2p

### ladder_2p/l1_n+8_l0_n+8

Parents: (1,+8) f 35.663 c/d, gamma -3.347e-10 s^-1; (0,+8) f 34.144 c/d, gamma -4.485e-07 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| harmonic | (1,+8) + (1,+8) -> (2,+17) | +1.731e-01 | +9.044e+01 |  | 5.396e-14 | +3.03e-05 | 1.432e+04 |
| sum | (0,+8) + (1,+8) -> (1,+17) | -1.098e-01 | +1.121e+02 |  | 2.217e-11 | +2.86e-05 | 1.890e+04 |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  | (1,+8) + (1,+8) -> (2,+17) | runaway | 2.57e+07 | 8.7255e-04 | 11.5 | 2.15e-12 | 1 |
| d2 | 2 | 2 |  | (1,+8) + (1,+8) -> (2,+17); (0,+8) + (1,+8) -> (1,+17) | runaway | 2.58e+07 | 8.6863e-04 | 11.5 | 6.973e-218 | 1 |

Same run from some N on: not reached.

### ladder_2p/l1_n+8_l1_n+7

Parents: (1,+8) f 35.663 c/d, gamma -3.347e-10 s^-1; (1,+7) f 32.096 c/d, gamma -7.808e-07 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| harmonic | (1,+8) + (1,+8) -> (2,+17) | +1.731e-01 | +9.044e+01 |  | 5.396e-14 | +3.03e-05 | 1.432e+04 |
| sum | (1,+7) + (1,+8) -> (2,+16) | -5.221e-02 | +1.006e+02 |  | 4.698e-11 | +2.67e-05 | 1.789e+04 |
| sum | (1,+7) + (1,+8) -> (0,+17) | +2.365e-01 | +1.033e+02 |  | 6.183e-11 | +2.69e-05 | 1.568e+04 |
| harmonic | (1,+7) + (1,+7) -> (0,+16) | +7.401e-03 | +1.147e+02 |  | 6.687e-11 | +2.30e-05 | 2.181e+04 |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  | (1,+8) + (1,+8) -> (2,+17) | runaway | 1.47e+07 | 8.7255e-04 | 11.5 | 1.503e-12 | 1 |
| d2 | 2 | 2 |  | (1,+8) + (1,+8) -> (2,+17); (1,+7) + (1,+8) -> (2,+16) | runaway | 1.48e+07 | 8.7327e-04 | 11.5 | 4.616e-125 | 1 |
| d3 | 3 | 3 |  | (1,+8) + (1,+8) -> (2,+17); (1,+7) + (1,+8) -> (2,+16); (1,+7) + (1,+8) -> (0,+17) | runaway | 1.48e+07 | 8.6697e-04 | 11.5 | 1.453e-130 | 2 |
| d4 | 4 | 4 |  | (1,+7) + (1,+7) -> (0,+16); (1,+8) + (1,+8) -> (2,+17); (1,+7) + (1,+8) -> (2,+16); (1,+7) + (1,+8) -> (0,+17) | end | 1e+13 | 6.6875e-11 | 5e-09 | 1 | 1 |

Same run from some N on: not reached.

### ladder_2p/l1_n+8_l2_n+7

Parents: (1,+8) f 35.663 c/d, gamma -3.347e-10 s^-1; (2,+7) f 33.770 c/d, gamma -5.357e-07 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| harmonic | (1,+8) + (1,+8) -> (2,+17) | +1.731e-01 | +9.044e+01 |  | 5.396e-14 | +3.03e-05 | 1.432e+04 |
| sum | (1,+8) + (2,+7) -> (3,+16) | -2.200e-01 | +1.029e+02 |  | 3.886e-11 | +2.82e-05 | 1.576e+04 |
| sum | (1,+8) + (2,+7) -> (1,+17) | +2.645e-01 | +9.027e+01 |  | 5.581e-11 | +2.85e-05 | 1.315e+04 |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  | (1,+8) + (1,+8) -> (2,+17) | runaway | 2.15e+07 | 8.7255e-04 | 11.5 | 1.936e-12 | 1 |
| d2 | 2 | 2 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,+7) -> (3,+16) | runaway | 2.16e+07 | 8.7209e-04 | 11.5 | 8.41e-173 | 1 |
| d3 | 3 | 3 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,+7) -> (3,+16); (1,+8) + (2,+7) -> (1,+17) | runaway | 2.16e+07 | 8.6930e-04 | 11.5 | 2.305e-204 | 2 |

Same run from some N on: not reached.

### ladder_2p/l2_n-2_l1_n+7

Parents: (1,+7) f 32.096 c/d, gamma -7.808e-07 s^-1; (2,-2) f 10.944 c/d, gamma -6.697e-11 s^-1

Built by: each parent's pairs ranked by E_th, the parents in turn; plus, for each combination of the parents, its direct triplet of highest R_c.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (2,-2) -> (3,-31) + (5,-8) | -4.038e-05 | +2.147e-01 | 4.870e-12 | 1.499e-10 | +1.47e-11 |  |
| pair | (2,-2) -> (5,-9) + (7,-50) | +8.218e-04 | +5.254e-02 | 1.302e-09 | 1.373e-09 | +2.48e-09 |  |
| pair | (2,-2) -> (4,-29) + (4,-7) | -1.662e-03 | +1.754e-01 | 8.448e-09 | 1.563e-07 | +2.03e-11 |  |
| pair | (2,-2) -> (4,-24) + (6,-11) | -1.126e-03 | +2.868e-01 | 1.088e-08 | 7.362e-09 | -3.67e-11 |  |
| pair | (2,-2) -> (5,-21) + (5,-11) | -8.884e-04 | +1.051e-01 | 7.174e-08 | 7.245e-09 | -5.08e-11 |  |
| pair | (2,-2) -> (6,-29) + (6,-12) | +2.325e-03 | +8.988e-02 | 1.762e-07 | 5.422e-06 | +1.47e-11 |  |
| pair | (2,-2) -> (6,-49) + (6,-10) | +7.677e-03 | +5.037e-02 | 2.559e-07 | 2.744e-07 | +1.89e-09 |  |
| pair | (2,-2) -> (6,-35) + (6,-11) | +3.389e-03 | +6.764e-02 | 3.602e-07 | 1.089e-06 | +9.07e-11 |  |
| pair | (2,-2) -> (6,-10) + (8,-64) | +1.326e-02 | +2.058e-02 | 7.872e-07 | 7.966e-07 | +1.13e-08 |  |
| harmonic | (1,+7) + (1,+7) -> (0,+16) | +7.401e-03 | +1.147e+02 |  |  | +2.30e-05 | (R_c 4.758e+08) |
| sum | (1,+7) + (2,-2) -> (1,+10) | -1.785e-02 | +4.909e+00 |  |  | +3.66e-06 | (R_c 4.404e+07) |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d4 | 4 | 3 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 1.9e-08 | 1 | 1 |
| d6 | 6 | 4 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 2.11e-08 | 1 | 1 |
| d8 | 8 | 5 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 2.14e-08 | 1 | 1 |
| d10 | 10 | 6 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 1.35e-08 | 1 | 1 |
| d12 | 12 | 7 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 2.03e-08 | 1 | 1 |
| d14 | 14 | 8 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 1.96e-08 | 1 | 1 |
| d16 | 16 | 9 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 1.63e-08 | 1 | 1 |
| d18 | 18 | 11 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 1.98e-08 | 1 | 1 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 4 on.

### ladder_2p/l2_n-2_l1_n+8

Parents: (2,-2) f 10.944 c/d, gamma -6.697e-11 s^-1; (1,+8) f 35.663 c/d, gamma -3.347e-10 s^-1

Built by: candidate triplets ranked by E_th (pairs) and E_eq (harmonic, sum), lowest first, each parent's lowest first.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| harmonic | (1,+8) + (1,+8) -> (2,+17) | +1.731e-01 | +9.044e+01 |  | 5.396e-14 | +3.03e-05 | 1.432e+04 |
| pair | (2,-2) -> (3,-31) + (5,-8) | -4.038e-05 | +2.147e-01 | 4.870e-12 | 1.499e-10 | +1.47e-11 |  |
| sum | (1,+8) + (2,-2) -> (1,+11) | +1.905e-01 | +4.973e+00 |  | 2.082e-11 | +7.59e-06 | 1.072e+03 |
| sum | (1,+8) + (2,-2) -> (3,+10) | -3.763e-01 | +5.533e+00 |  | 5.809e-11 | +7.13e-06 | 6.578e+02 |
| pair | (2,-2) -> (5,-9) + (7,-50) | +8.218e-04 | +5.254e-02 | 1.302e-09 | 1.373e-09 | +2.48e-09 |  |
| pair | (2,-2) -> (4,-29) + (4,-7) | -1.662e-03 | +1.754e-01 | 8.448e-09 | 1.563e-07 | +2.03e-11 |  |
| pair | (2,-2) -> (4,-24) + (6,-11) | -1.126e-03 | +2.868e-01 | 1.088e-08 | 7.362e-09 | -3.67e-11 |  |
| pair | (2,-2) -> (5,-21) + (5,-11) | -8.884e-04 | +1.051e-01 | 7.174e-08 | 7.245e-09 | -5.08e-11 |  |
| pair | (2,-2) -> (6,-29) + (6,-12) | +2.325e-03 | +8.988e-02 | 1.762e-07 | 5.422e-06 | +1.47e-11 |  |
| pair | (2,-2) -> (6,-49) + (6,-10) | +7.677e-03 | +5.037e-02 | 2.559e-07 | 2.744e-07 | +1.89e-09 |  |
| pair | (2,-2) -> (6,-35) + (6,-11) | +3.389e-03 | +6.764e-02 | 3.602e-07 | 1.089e-06 | +9.07e-11 |  |
| pair | (2,-2) -> (3,-6) + (5,-29) | +1.660e-02 | +1.748e-01 | 6.420e-07 | 1.674e-05 | +1.63e-11 |  |
| pair | (2,-2) -> (4,-13) + (4,-11) | -2.590e-03 | +1.458e-01 | 6.569e-07 | 2.949e-09 | -6.28e-11 |  |
| pair | (2,-2) -> (6,-10) + (8,-64) | +1.326e-02 | +2.058e-02 | 7.872e-07 | 7.966e-07 | +1.13e-08 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 1 | 1 |  | (1,+8) + (1,+8) -> (2,+17) | runaway | 1.72e+11 | 8.7486e-04 | 11.5 | 6.878e-11 | 1 |
| d3 | 3 | 2 |  | (1,+8) + (1,+8) -> (2,+17) | end | 1e+13 | 1.4984e-10 | 0.00843 | 1 | 2 |
| d4 | 4 | 3 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,-2) -> (1,+11) | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d5 | 5 | 4 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,-2) -> (1,+11); (1,+8) + (2,-2) -> (3,+10) | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d7 | 7 | 5 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,-2) -> (1,+11); (1,+8) + (2,-2) -> (3,+10) | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d9 | 9 | 6 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,-2) -> (1,+11); (1,+8) + (2,-2) -> (3,+10) | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d11 | 11 | 7 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,-2) -> (1,+11); (1,+8) + (2,-2) -> (3,+10) | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d13 | 13 | 8 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,-2) -> (1,+11); (1,+8) + (2,-2) -> (3,+10) | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d15 | 15 | 9 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,-2) -> (1,+11); (1,+8) + (2,-2) -> (3,+10) | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d17 | 17 | 10 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,-2) -> (1,+11); (1,+8) + (2,-2) -> (3,+10) | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d18 | 18 | 11 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,-2) -> (1,+11); (1,+8) + (2,-2) -> (3,+10) | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d20 | 20 | 12 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,-2) -> (1,+11); (1,+8) + (2,-2) -> (3,+10) | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d22 | 22 | 13 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,-2) -> (1,+11); (1,+8) + (2,-2) -> (3,+10) | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |
| d23 | 23 | 14 |  | (1,+8) + (1,+8) -> (2,+17); (1,+8) + (2,-2) -> (1,+11); (1,+8) + (2,-2) -> (3,+10) | end | 1e+13 | 1.4979e-10 | 0.00843 | 1 | 2 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 3 on.

### ladder_2p/l2_n-3_l2_n-2

Parents: (2,-3) f 9.373 c/d, gamma -8.642e-13 s^-1; (2,-2) f 10.944 c/d, gamma -6.697e-11 s^-1

Built by: step k: the k daughters shared by both parents with the lowest larger E_th of their two pairs, each with both pairs.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | gamma sum [1/s] | mu |
|---|---|---|---|---|---|---|---|
| pair | (2,-3) -> (6,-55) + (6,-12) | +9.957e-04 | +4.810e-03 | 3.607e-07 | 3.608e-07 | +4.83e-09 |  |
| pair | (2,-2) -> (6,-29) + (6,-12) | +2.325e-03 | +8.988e-02 | 1.762e-07 | 5.422e-06 | +1.47e-11 |  |
| pair | (2,-3) -> (1,-16) + (3,-6) | -4.566e-03 | +1.934e-01 | 1.328e-06 | 1.973e-06 | +3.94e-12 |  |
| pair | (2,-2) -> (3,-6) + (5,-29) | +1.660e-02 | +1.748e-01 | 6.420e-07 | 1.674e-05 | +1.63e-11 |  |
| pair | (2,-3) -> (3,-12) + (3,-10) | +1.305e-03 | -2.282e-01 | 8.862e-08 | 1.552e-07 | +2.67e-12 |  |
| pair | (2,-2) -> (3,-10) + (5,-14) | -8.906e-03 | +1.300e-01 | 7.776e-06 | 5.008e-08 | -6.20e-11 |  |
| pair | (2,-3) -> (6,-13) + (8,-62) | -8.013e-03 | -7.567e-03 | 4.409e-06 | 4.410e-06 | +1.08e-08 |  |
| pair | (2,-2) -> (4,-18) + (6,-13) | +1.936e-02 | +3.112e-01 | 7.792e-06 | 1.549e-07 | -5.87e-11 |  |
| pair | (2,-3) -> (6,-17) + (8,-35) | +9.380e-04 | -1.569e-01 | 1.658e-08 | 1.680e-08 | +1.35e-10 |  |
| pair | (2,-2) -> (6,-13) + (8,-35) | -1.191e-02 | +5.356e-02 | 9.060e-06 | 3.676e-05 | +6.61e-11 |  |
| pair | (2,-3) -> (3,-23) + (5,-11) | -1.354e-02 | +1.473e-01 | 1.175e-05 | 1.283e-05 | +1.92e-11 |  |
| pair | (2,-2) -> (3,-23) + (5,-9) | -2.944e-02 | +3.069e-01 | 6.225e-06 | 9.885e-07 | -4.79e-11 |  |
| pair | (2,-2) -> (5,-21) + (5,-11) | -8.884e-04 | +1.051e-01 | 7.174e-08 | 7.245e-09 | -5.08e-11 |  |
| pair | (2,-3) -> (4,-26) + (6,-14) | +7.398e-04 | +1.153e-01 | 4.925e-08 | 5.140e-08 | +4.07e-11 |  |
| pair | (2,-2) -> (6,-14) + (8,-31) | -7.186e-03 | +4.374e-02 | 1.192e-05 | 1.216e-03 | +7.36e-12 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | <D>_T/<P>_T | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d3 | 3 | 2 | (6,-12) |  | end | 1e+13 | 5.3749e-06 | 5.31e-08 | 0.9997 | 2 |
| d6 | 6 | 4 | (3,-6) (6,-12) |  | end | 1e+13 | 1.2484e-05 | 4.3 | 0.7838 | 4 |
| d9 | 9 | 6 | (3,-6) (3,-10) (6,-12) |  | runaway | 1.06e+12 | 1.9509e-04 | 9.7 | 0.2884 | 5 |
| d12 | 12 | 8 | (3,-6) (3,-10) (6,-12) (6,-13) |  | runaway | 1.01e+12 | 1.6660e-04 | 20.8 | 0.1023 | 5 |
| d14 | 14 | 10 | (3,-6) (3,-10) (6,-12) (8,-35) (6,-13) |  | runaway | 1.5e+12 | 1.2066e-04 | 26.5 | 0.1424 | 6 |
| d17 | 17 | 12 | (3,-6) (3,-10) (3,-23) (6,-12) (8,-35) (6,-13) |  | runaway | 1.5e+12 | 1.1786e-04 | 29.7 | 0.1477 | 6 |
| d18 | 18 | 13 | (3,-6) (3,-10) (3,-23) (5,-11) (6,-12) (8,-35) (6,-13) |  | end | 1e+13 | 1.4526e-06 | 4.73 | 0.7866 | 9 |
| d21 | 21 | 15 | (3,-6) (3,-10) (3,-23) (5,-11) (6,-14) (6,-12) (8,-35) (6,-13) |  | end | 1e+13 | 1.4110e-06 | 2.4 | 0.9987 | 17 |

Same run from some N on: not reached.

