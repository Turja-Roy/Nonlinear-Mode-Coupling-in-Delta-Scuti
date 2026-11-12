# dsct_M1.85_T7555: ladders and runs

Written by scripts/network_index.py.

Notation. (l,n): degree and radial order of a mode; f: frequency; gamma: linear rate, < 0 driven (parent),
> 0 damped (daughter). N: number of daughters. Pair p -> b + c: parametric triplet; harmonic x + x -> c and
sum x + y -> c: direct triplets. Delta: detuning; kappa: coupling coefficient. E_th: parent energy above
which a pair grows (Weinberg+12 eq. 50); E_eq: parent energy at a triplet's fixed point (MW25 eq. 6), for a
sum triplet E_x + E_y. Energies in E_star = GM^2/R. Shared daughter: coupled to both parents. Run columns,
second half of the run: <E_par>, mean summed parent energy; swing, (max - min)/mean of it;
dissipation/driving, sum 2 gamma_d <E_d> over sum 2|gamma_p| <E_p>; n_active, daughters above 1e-3 of the
largest. stop: end (t = 1e13 s), runaway (summed energy 1e-2), wall (time limit).

## ladder_1p

### ladder_1p/l2_n-4

Parents: (2,-4) f 8.193 c/d, gamma -1.883e-14 s^-1

Built by: pairs ranked by E_th, lowest first; a pair that would bring N past the target skipped.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (2,-4) -> (8,-32) + (10,-43) | +1.232e-04 | -1.494e-01 | 6.585e-10 | 6.585e-10 |  |
| pair | (2,-4) -> (11,-48) + (13,-51) | -3.369e-04 | +2.375e-01 | 8.515e-09 | 8.515e-09 |  |
| pair | (2,-4) -> (13,-55) + (15,-60) | -9.916e-04 | -6.227e-01 | 9.445e-09 | 9.445e-09 |  |
| pair | (2,-4) -> (10,-52) + (12,-41) | -1.181e-04 | +6.082e-02 | 1.607e-08 | 1.607e-08 |  |
| pair | (2,-4) -> (9,-45) + (11,-39) | -8.001e-04 | -2.005e-01 | 4.749e-08 | 4.749e-08 |  |
| pair | (2,-4) -> (2,-14) + (2,-7) | -2.161e-03 | +2.574e-01 | 1.475e-07 | 1.494e-07 |  |
| pair | (2,-4) -> (12,-53) + (14,-54) | +1.164e-03 | +1.863e-01 | 1.514e-07 | 1.514e-07 |  |
| pair | (2,-4) -> (8,-60) + (10,-28) | +7.441e-04 | +1.799e-02 | 1.730e-07 | 1.730e-07 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | end | 1e+13 | 1.3283e-12 | 0.188 | 0.0000 | 0 |
| d4 | 4 | 2 |  |  | end | 1e+13 | 1.3283e-12 | 0.188 | 0.0000 | 0 |
| d6 | 6 | 3 |  |  | end | 1e+13 | 1.3283e-12 | 0.188 | 0.0000 | 0 |
| d8 | 8 | 4 |  |  | end | 1e+13 | 1.3283e-12 | 0.188 | 0.0000 | 0 |
| d10 | 10 | 5 |  |  | end | 1e+13 | 1.3283e-12 | 0.188 | 0.0000 | 0 |
| d12 | 12 | 6 |  |  | end | 1e+13 | 1.3283e-12 | 0.188 | 0.0000 | 1 |
| d14 | 14 | 7 |  |  | end | 1e+13 | 1.3283e-12 | 0.188 | 0.0000 | 1 |
| d16 | 16 | 8 |  |  | end | 1e+13 | 1.3283e-12 | 0.188 | 0.0000 | 1 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 12 on.

## ladder_2p

### ladder_2p/l0_n+1_l2_n+6

Parents: (0,+1) f 12.557 c/d, gamma -1.702e-10 s^-1; (2,+6) f 35.374 c/d, gamma -1.903e-07 s^-1

Built by: each parent's pairs ranked by E_th, the parents in turn; plus, for each combination of the parents, its direct triplet of highest R_c.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (0,+1) -> (5,-31) + (5,-8) | +7.693e-03 | +1.164e-01 | 1.277e-07 | none |  |
| pair | (0,+1) -> (4,-22) + (4,-7) | +1.750e-02 | +2.558e-01 | 4.258e-07 | none |  |
| pair | (0,+1) -> (3,-12) + (3,-6) | -5.042e-03 | +1.933e-01 | 5.151e-07 | none |  |
| pair | (0,+1) -> (6,-41) + (6,-9) | -3.190e-02 | +6.725e-02 | 1.350e-06 | 1.663e-06 |  |
| pair | (0,+1) -> (6,-24) + (6,-11) | -1.732e-02 | +1.955e-01 | 3.221e-06 | none |  |
| pair | (0,+1) -> (5,-18) + (5,-10) | -1.376e-02 | +2.257e-01 | 3.278e-06 | none |  |
| pair | (0,+1) -> (6,-40) + (6,-9) | +3.679e-02 | +5.369e-02 | 3.488e-06 | 4.567e-06 |  |
| pair | (0,+1) -> (4,-16) + (4,-8) | +2.248e-02 | +2.527e-01 | 5.604e-06 | none |  |
| pair | (0,+1) -> (5,-32) + (5,-8) | -6.548e-02 | +1.264e-01 | 8.658e-06 | none |  |
| sum | (0,+1) + (2,+6) -> (2,+9) | +5.251e-02 | +5.332e+00 |  |  | (R_c 4.492e+07) |
| harmonic | (2,+6) + (2,+6) -> (0,+15) | -4.083e-02 | +1.217e+02 |  |  | (R_c 7.649e+08) |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d4 | 4 | 3 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 1.1088e-11 | 1.25e-08 | 0.9994 | 1 |
| d6 | 6 | 4 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 1.1088e-11 | 1.47e-09 | 0.9994 | 1 |
| d8 | 8 | 5 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 1.1088e-11 | 1.31e-12 | 0.9994 | 1 |
| d10 | 10 | 6 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 1.1088e-11 | 2.28e-12 | 0.9994 | 1 |
| d12 | 12 | 7 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 1.1088e-11 | 1.06e-12 | 0.9994 | 1 |
| d14 | 14 | 8 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 1.1088e-11 | 3.64e-12 | 0.9994 | 1 |
| d16 | 16 | 10 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 1.1088e-11 | 1.14e-09 | 0.9994 | 1 |
| d18 | 18 | 11 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 1.1088e-11 | 2.98e-12 | 0.9994 | 1 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 4 on.

### ladder_2p/l1_n-1_l2_n+4

Parents: (1,-1) f 12.839 c/d, gamma -2.211e-10 s^-1; (2,+4) f 28.090 c/d, gamma -3.122e-08 s^-1

Built by: each parent's pairs ranked by E_th, the parents in turn; plus, for each combination of the parents, its direct triplet of highest R_c.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (1,-1) -> (6,-9) + (7,-42) | +1.603e-03 | +8.654e-02 | 1.580e-09 | 1.997e-09 |  |
| pair | (1,-1) -> (5,-8) + (6,-34) | -1.750e-02 | +1.596e-01 | 1.473e-07 | 1.084e-06 |  |
| pair | (1,-1) -> (4,-23) + (5,-8) | -1.345e-02 | +4.045e-01 | 1.947e-07 | none |  |
| pair | (1,-1) -> (4,-14) + (5,-10) | +2.454e-03 | +1.994e-01 | 2.115e-07 | none |  |
| pair | (1,-1) -> (5,-30) + (6,-9) | +3.584e-02 | +2.498e-01 | 1.397e-06 | none |  |
| pair | (1,-1) -> (4,-7) + (5,-25) | -4.419e-02 | +2.276e-01 | 1.699e-06 | none |  |
| pair | (1,-1) -> (6,-12) + (7,-23) | +8.770e-03 | +1.338e-01 | 1.729e-06 | none |  |
| pair | (1,-1) -> (6,-9) + (7,-43) | -6.671e-02 | +7.783e-02 | 3.219e-06 | 4.001e-06 |  |
| pair | (1,-1) -> (6,-9) + (7,-41) | +7.688e-02 | +9.461e-02 | 3.554e-06 | 4.727e-06 |  |
| pair | (1,-1) -> (5,-31) + (6,-9) | -6.016e-02 | +2.259e-01 | 4.133e-06 | none |  |
| pair | (1,-1) -> (5,-8) + (6,-33) | +7.381e-02 | +1.535e-01 | 4.295e-06 | 2.382e-02 |  |
| sum | (1,-1) + (2,+4) -> (3,+7) | +7.466e-02 | +1.375e+00 |  |  | (R_c 2.229e+06) |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d3 | 3 | 2 |  | (1,-1) + (2,+4) -> (3,+7) | runaway | 3.69e+08 | 8.7214e-04 | 11.5 | 0.0000 | 3 |
| d5 | 5 | 3 |  | (1,-1) + (2,+4) -> (3,+7) | runaway | 3.69e+08 | 8.7211e-04 | 11.5 | 0.0000 | 5 |
| d7 | 7 | 5 |  | (1,-1) + (2,+4) -> (3,+7) | runaway | 3.69e+08 | 8.7222e-04 | 11.5 | 0.0000 | 7 |
| d9 | 9 | 6 |  | (1,-1) + (2,+4) -> (3,+7) | runaway | 3.69e+08 | 8.7212e-04 | 11.5 | 0.0000 | 9 |
| d11 | 11 | 7 |  | (1,-1) + (2,+4) -> (3,+7) | runaway | 3.69e+08 | 8.7211e-04 | 11.5 | 0.0000 | 11 |
| d13 | 13 | 8 |  | (1,-1) + (2,+4) -> (3,+7) | runaway | 3.69e+08 | 8.7213e-04 | 11.5 | 0.0000 | 13 |
| d15 | 15 | 10 |  | (1,-1) + (2,+4) -> (3,+7) | runaway | 3.69e+08 | 8.7211e-04 | 11.5 | 0.0000 | 15 |
| d17 | 17 | 12 |  | (1,-1) + (2,+4) -> (3,+7) | runaway | 3.69e+08 | 8.7222e-04 | 11.5 | 0.0000 | 17 |

Same run from some N on: not reached.

### ladder_2p/l2_n-4_l1_n-1

Parents: (1,-1) f 12.839 c/d, gamma -2.211e-10 s^-1; (2,-4) f 8.193 c/d, gamma -1.883e-14 s^-1

Built by: step k: the k daughters shared by both parents with the lowest larger E_th of their two pairs, each with both pairs.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (1,-1) -> (5,-8) + (6,-33) | +7.381e-02 | +1.535e-01 | 4.295e-06 | 2.382e-02 |  |
| pair | (2,-4) -> (6,-33) + (8,-27) | +3.625e-03 | -1.839e-01 | 1.016e-06 | 1.016e-06 |  |
| pair | (1,-1) -> (6,-10) + (7,-31) | +3.558e-02 | +1.348e-01 | 5.169e-06 | none |  |
| pair | (2,-4) -> (5,-20) + (7,-31) | +1.117e-02 | -1.523e-01 | 4.360e-06 | 4.361e-06 |  |
| pair | (1,-1) -> (2,-7) + (3,-7) | +3.141e-03 | +5.148e-02 | 5.789e-06 | none |  |
| pair | (2,-4) -> (2,-14) + (2,-7) | -2.161e-03 | +2.574e-01 | 1.475e-07 | 1.494e-07 |  |
| pair | (1,-1) -> (6,-9) + (7,-43) | -6.671e-02 | +7.783e-02 | 3.219e-06 | 4.001e-06 |  |
| pair | (2,-4) -> (7,-43) + (9,-28) | -6.193e-03 | +4.745e-02 | 1.036e-05 | 1.036e-05 |  |
| pair | (1,-1) -> (6,-9) + (7,-42) | +1.603e-03 | +8.654e-02 | 1.580e-09 | 1.997e-09 |  |
| pair | (2,-4) -> (5,-16) + (7,-42) | +2.011e-03 | +4.186e-03 | 1.112e-05 | 1.112e-05 |  |
| pair | (1,-1) -> (5,-8) + (6,-36) | -1.994e-01 | +1.450e-01 | 1.359e-05 | 3.201e-05 |  |
| pair | (2,-4) -> (6,-36) + (6,-19) | +2.039e-03 | +1.047e-02 | 7.851e-06 | 7.852e-06 |  |
| pair | (1,-1) -> (4,-7) + (5,-24) | +9.638e-02 | +2.659e-01 | 6.358e-06 | none |  |
| pair | (2,-4) -> (5,-24) + (7,-26) | +7.906e-03 | -1.223e-01 | 1.570e-05 | 1.571e-05 |  |
| pair | (1,-1) -> (6,-11) + (7,-26) | +1.468e-02 | +1.351e-01 | 4.294e-06 | none |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d3 | 3 | 2 | (6,-33) |  | runaway | 6.37e+10 | 8.4818e-04 | 2.77 | 0.1663 | 2 |
| d6 | 6 | 4 | (7,-31) (6,-33) |  | runaway | 7.2e+10 | 7.8657e-04 | 3.05 | 0.2034 | 2 |
| d9 | 9 | 6 | (2,-7) (7,-31) (6,-33) |  | runaway | 7.99e+10 | 7.0479e-04 | 3.37 | 0.2011 | 5 |
| d12 | 12 | 8 | (2,-7) (7,-31) (6,-33) (7,-43) |  | end | 1e+13 | 3.7025e-06 | 0.00338 | 1.0009 | 6 |
| d14 | 14 | 10 | (2,-7) (7,-42) (7,-31) (6,-33) (7,-43) |  | end | 1e+13 | 2.0088e-09 | 3.9 | 0.9911 | 1 |
| d16 | 16 | 12 | (2,-7) (7,-42) (7,-31) (6,-36) (6,-33) (7,-43) |  | end | 1e+13 | 2.0118e-09 | 3.88 | 0.9902 | 1 |
| d19 | 19 | 14 | (2,-7) (7,-42) (5,-24) (7,-31) (6,-36) (6,-33) (7,-43) |  | end | 1e+13 | 2.0149e-09 | 3.87 | 0.9893 | 1 |
| d20 | 20 | 15 | (2,-7) (7,-26) (7,-42) (5,-24) (7,-31) (6,-36) (6,-33) (7,-43) |  | end | 1e+13 | 2.0145e-09 | 3.87 | 0.9894 | 1 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 19 on.

