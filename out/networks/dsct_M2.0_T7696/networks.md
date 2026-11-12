# dsct_M2.0_T7696: ladders and runs

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

### ladder_1p/l2_n-2

Parents: (2,-2) f 10.944 c/d, gamma -6.697e-11 s^-1

Built by: pairs ranked by E_th, lowest first; a pair that would bring N past the target skipped.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (2,-2) -> (3,-31) + (5,-8) | -4.038e-05 | +2.147e-01 | 4.870e-12 | 1.499e-10 |  |
| pair | (2,-2) -> (5,-9) + (7,-50) | +8.218e-04 | +5.254e-02 | 1.302e-09 | 1.373e-09 |  |
| pair | (2,-2) -> (4,-29) + (4,-7) | -1.662e-03 | +1.754e-01 | 8.448e-09 | 1.563e-07 |  |
| pair | (2,-2) -> (4,-24) + (6,-11) | -1.126e-03 | +2.868e-01 | 1.088e-08 | none |  |
| pair | (2,-2) -> (5,-21) + (5,-11) | -8.884e-04 | +1.051e-01 | 7.174e-08 | none |  |
| pair | (2,-2) -> (6,-29) + (6,-12) | +2.325e-03 | +8.988e-02 | 1.762e-07 | 5.421e-06 |  |
| pair | (2,-2) -> (6,-49) + (6,-10) | +7.677e-03 | +5.037e-02 | 2.559e-07 | 2.744e-07 |  |
| pair | (2,-2) -> (6,-35) + (6,-11) | +3.389e-03 | +6.764e-02 | 3.602e-07 | 1.089e-06 |  |
| pair | (2,-2) -> (6,-10) + (8,-64) | +1.326e-02 | +2.058e-02 | 7.872e-07 | 7.965e-07 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | end | 1e+13 | 1.4977e-10 | 0.00843 | 1.0005 | 2 |
| d4 | 4 | 2 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1.0005 | 2 |
| d6 | 6 | 3 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1.0005 | 2 |
| d8 | 8 | 4 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1.0005 | 2 |
| d10 | 10 | 5 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1.0005 | 2 |
| d12 | 12 | 6 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1.0005 | 2 |
| d14 | 14 | 7 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1.0005 | 2 |
| d16 | 16 | 9 |  |  | end | 1e+13 | 1.4979e-10 | 0.00843 | 1.0005 | 2 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 2 on.

## ladder_2p

### ladder_2p/l2_n-2_l1_n+7

Parents: (1,+7) f 32.096 c/d, gamma -7.808e-07 s^-1; (2,-2) f 10.944 c/d, gamma -6.697e-11 s^-1

Built by: each parent's pairs ranked by E_th, the parents in turn; plus, for each combination of the parents, its direct triplet of highest R_c.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (2,-2) -> (3,-31) + (5,-8) | -4.038e-05 | +2.147e-01 | 4.870e-12 | 1.499e-10 |  |
| pair | (2,-2) -> (5,-9) + (7,-50) | +8.218e-04 | +5.254e-02 | 1.302e-09 | 1.373e-09 |  |
| pair | (2,-2) -> (4,-29) + (4,-7) | -1.662e-03 | +1.754e-01 | 8.448e-09 | 1.563e-07 |  |
| pair | (2,-2) -> (4,-24) + (6,-11) | -1.126e-03 | +2.868e-01 | 1.088e-08 | none |  |
| pair | (2,-2) -> (5,-21) + (5,-11) | -8.884e-04 | +1.051e-01 | 7.174e-08 | none |  |
| pair | (2,-2) -> (6,-29) + (6,-12) | +2.325e-03 | +8.988e-02 | 1.762e-07 | 5.421e-06 |  |
| pair | (2,-2) -> (6,-49) + (6,-10) | +7.677e-03 | +5.037e-02 | 2.559e-07 | 2.744e-07 |  |
| pair | (2,-2) -> (6,-35) + (6,-11) | +3.389e-03 | +6.764e-02 | 3.602e-07 | 1.089e-06 |  |
| pair | (2,-2) -> (6,-10) + (8,-64) | +1.326e-02 | +2.058e-02 | 7.872e-07 | 7.965e-07 |  |
| harmonic | (1,+7) + (1,+7) -> (0,+16) | +7.401e-03 | +1.147e+02 |  |  | (R_c 4.758e+08) |
| sum | (1,+7) + (2,-2) -> (1,+10) | -1.785e-02 | +4.909e+00 |  |  | (R_c 4.404e+07) |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d4 | 4 | 3 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 1.9e-08 | 1.0001 | 1 |
| d6 | 6 | 4 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 2.11e-08 | 1.0001 | 1 |
| d8 | 8 | 5 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 2.14e-08 | 1.0001 | 1 |
| d10 | 10 | 6 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 1.35e-08 | 1.0001 | 1 |
| d12 | 12 | 7 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 2.03e-08 | 1.0001 | 1 |
| d14 | 14 | 8 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 1.96e-08 | 1.0001 | 1 |
| d16 | 16 | 9 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 1.63e-08 | 1.0001 | 1 |
| d18 | 18 | 11 |  | (1,+7) + (1,+7) -> (0,+16); (1,+7) + (2,-2) -> (1,+10) | end | 1e+13 | 6.6875e-11 | 1.98e-08 | 1.0001 | 1 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 4 on.

### ladder_2p/l2_n-3_l2_n-2

Parents: (2,-3) f 9.373 c/d, gamma -8.642e-13 s^-1; (2,-2) f 10.944 c/d, gamma -6.697e-11 s^-1

Built by: step k: the k daughters shared by both parents with the lowest larger E_th of their two pairs, each with both pairs.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (2,-3) -> (6,-55) + (6,-12) | +9.957e-04 | +4.810e-03 | 3.607e-07 | 3.608e-07 |  |
| pair | (2,-2) -> (6,-29) + (6,-12) | +2.325e-03 | +8.988e-02 | 1.762e-07 | 5.421e-06 |  |
| pair | (2,-3) -> (1,-16) + (3,-6) | -4.566e-03 | +1.934e-01 | 1.328e-06 | 1.974e-06 |  |
| pair | (2,-2) -> (3,-6) + (5,-29) | +1.660e-02 | +1.748e-01 | 6.420e-07 | 1.673e-05 |  |
| pair | (2,-3) -> (3,-12) + (3,-10) | +1.305e-03 | -2.282e-01 | 8.862e-08 | 1.552e-07 |  |
| pair | (2,-2) -> (3,-10) + (5,-14) | -8.906e-03 | +1.300e-01 | 7.776e-06 | none |  |
| pair | (2,-3) -> (6,-13) + (8,-62) | -8.013e-03 | -7.567e-03 | 4.409e-06 | 4.410e-06 |  |
| pair | (2,-2) -> (4,-18) + (6,-13) | +1.936e-02 | +3.112e-01 | 7.792e-06 | none |  |
| pair | (2,-3) -> (6,-17) + (8,-35) | +9.380e-04 | -1.569e-01 | 1.658e-08 | 1.680e-08 |  |
| pair | (2,-2) -> (6,-13) + (8,-35) | -1.191e-02 | +5.356e-02 | 9.060e-06 | 3.674e-05 |  |
| pair | (2,-3) -> (3,-23) + (5,-11) | -1.354e-02 | +1.473e-01 | 1.175e-05 | 1.283e-05 |  |
| pair | (2,-2) -> (3,-23) + (5,-9) | -2.944e-02 | +3.069e-01 | 6.225e-06 | none |  |
| pair | (2,-2) -> (5,-21) + (5,-11) | -8.884e-04 | +1.051e-01 | 7.174e-08 | none |  |
| pair | (2,-3) -> (4,-26) + (6,-14) | +7.398e-04 | +1.153e-01 | 4.925e-08 | 5.137e-08 |  |
| pair | (2,-2) -> (6,-14) + (8,-31) | -7.186e-03 | +4.374e-02 | 1.192e-05 | 1.216e-03 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d3 | 3 | 2 | (6,-12) |  | end | 1e+13 | 5.3749e-06 | 5.31e-08 | 0.9997 | 2 |
| d6 | 6 | 4 | (3,-6) (6,-12) |  | end | 1e+13 | 1.2484e-05 | 4.3 | 0.8186 | 4 |
| d9 | 9 | 6 | (3,-6) (3,-10) (6,-12) |  | runaway | 1.06e+12 | 1.9509e-04 | 9.7 | 0.2831 | 5 |
| d12 | 12 | 8 | (3,-6) (3,-10) (6,-12) (6,-13) |  | runaway | 1.01e+12 | 1.6660e-04 | 20.8 | 0.1067 | 5 |
| d14 | 14 | 10 | (3,-6) (3,-10) (6,-12) (8,-35) (6,-13) |  | runaway | 1.5e+12 | 1.2066e-04 | 26.5 | 0.1490 | 6 |
| d17 | 17 | 12 | (3,-6) (3,-10) (3,-23) (6,-12) (8,-35) (6,-13) |  | runaway | 1.5e+12 | 1.1786e-04 | 29.7 | 0.1540 | 6 |
| d18 | 18 | 13 | (3,-6) (3,-10) (3,-23) (5,-11) (6,-12) (8,-35) (6,-13) |  | end | 1e+13 | 1.4526e-06 | 4.73 | 0.9595 | 9 |
| d21 | 21 | 15 | (3,-6) (3,-10) (3,-23) (5,-11) (6,-14) (6,-12) (8,-35) (6,-13) |  | end | 1e+13 | 1.4110e-06 | 2.4 | 0.9858 | 17 |

Same run from some N on: not reached.

