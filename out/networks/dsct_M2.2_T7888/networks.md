# dsct_M2.2_T7888: ladders and runs

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

### ladder_1p/l1_n+3

Parents: (1,+3) f 15.974 c/d, gamma -5.538e-09 s^-1

Built by: pairs ranked by E_th, lowest first; a pair that would bring N past the target skipped.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (1,+3) -> (1,-4) + (2,-1) | +2.360e-03 | +1.500e-01 | 3.593e-10 | none |  |
| pair | (1,+3) -> (5,-2) + (6,-67) | -2.503e-03 | +1.797e+00 | 3.653e-10 | 3.035e-07 |  |
| pair | (1,+3) -> (6,-5) + (7,-29) | -6.572e-03 | +1.945e+00 | 4.755e-10 | none |  |
| pair | (1,+3) -> (5,-26) + (6,-4) | -3.522e-04 | +2.113e-01 | 2.999e-09 | none |  |
| pair | (1,+3) -> (1,-14) + (2,+1) | +2.560e-03 | +4.848e-01 | 3.021e-09 | none |  |
| pair | (1,+3) -> (5,-21) + (6,-5) | -2.094e-02 | +1.257e+00 | 3.821e-09 | none |  |
| pair | (1,+3) -> (4,-46) + (5,-2) | +2.679e-03 | +1.443e+00 | 4.045e-09 | none |  |
| pair | (1,+3) -> (3,+0) + (4,-41) | -4.766e-03 | +2.332e+00 | 9.719e-09 | none |  |
| pair | (1,+3) -> (4,-1) + (5,-58) | +1.163e-02 | +4.294e+00 | 1.483e-08 | none |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | runaway | 2.97e+09 | 6.9850e-04 | 13.1 | 0.1041 | 2 |
| d4 | 4 | 2 |  |  | end | 1e+13 | 3.0096e-07 | 4.02e-09 | 0.9998 | 2 |
| d6 | 6 | 3 |  |  | runaway | 1.83e+10 | 2.7713e-04 | 16.5 | 0.7221 | 3 |
| d8 | 8 | 4 |  |  | runaway | 1.83e+10 | 2.8571e-04 | 15.5 | 0.6795 | 3 |
| d10 | 10 | 5 |  |  | runaway | 1.94e+10 | 2.8613e-04 | 15 | 0.6356 | 3 |
| d12 | 12 | 7 |  |  | runaway | 1.92e+10 | 3.0274e-04 | 11.5 | 0.7261 | 4 |
| d14 | 14 | 8 |  |  | runaway | 1.7e+10 | 3.1058e-04 | 13.6 | 0.6786 | 4 |
| d16 | 16 | 9 |  |  | end | 1e+13 | 2.2734e-09 | 1.33 | 0.9996 | 4 |

Same run from some N on: not reached.

## ladder_2p

### ladder_2p/l1_n+3_l1_n+4

Parents: (1,+3) f 15.974 c/d, gamma -5.538e-09 s^-1; (1,+4) f 17.822 c/d, gamma -7.208e-08 s^-1

Built by: each parent's pairs ranked by E_th, the parents in turn; plus, for each combination of the parents, its direct triplet of highest R_c.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (1,+3) -> (1,-4) + (2,-1) | +2.360e-03 | +1.500e-01 | 3.593e-10 | none |  |
| pair | (1,+4) -> (0,+1) + (1,-1) | +1.092e-03 | +8.791e-01 | 1.196e-09 | none |  |
| pair | (1,+3) -> (5,-2) + (6,-67) | -2.503e-03 | +1.797e+00 | 3.653e-10 | 3.035e-07 |  |
| pair | (1,+4) -> (3,-15) + (4,-1) | +3.022e-02 | +9.673e-01 | 3.703e-08 | none |  |
| pair | (1,+3) -> (6,-5) + (7,-29) | -6.572e-03 | +1.945e+00 | 4.755e-10 | none |  |
| pair | (1,+4) -> (2,+1) + (3,-15) | -1.627e-02 | +6.737e-01 | 3.939e-08 | none |  |
| pair | (1,+3) -> (5,-26) + (6,-4) | -3.522e-04 | +2.113e-01 | 2.999e-09 | none |  |
| pair | (1,+4) -> (4,-1) + (5,-25) | -6.620e-02 | +2.824e+00 | 1.149e-07 | none |  |
| pair | (1,+3) -> (1,-14) + (2,+1) | +2.560e-03 | +4.848e-01 | 3.021e-09 | none |  |
| pair | (1,+4) -> (4,-1) + (5,-24) | +8.358e-02 | +2.741e+00 | 1.659e-07 | none |  |
| sum | (1,+3) + (1,+4) -> (2,+9) | +1.995e-03 | +8.019e+00 |  |  | (R_c 6.170e+07) |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d5 | 5 | 3 |  | (1,+3) + (1,+4) -> (2,+9) | runaway | 1.85e+08 | 7.6372e-04 | 13.1 | 0.0022 | 2 |
| d7 | 7 | 4 |  | (1,+3) + (1,+4) -> (2,+9) | runaway | 1.85e+08 | 7.6259e-04 | 13.1 | 0.0023 | 2 |
| d9 | 9 | 5 |  | (1,+3) + (1,+4) -> (2,+9) | runaway | 2.15e+08 | 4.5896e-04 | 15 | 0.0051 | 4 |
| d11 | 11 | 6 |  | (1,+3) + (1,+4) -> (2,+9) | runaway | 2.1e+08 | 6.6907e-04 | 14.7 | 0.0041 | 4 |
| d13 | 13 | 8 |  | (1,+3) + (1,+4) -> (2,+9) | runaway | 2.52e+08 | 5.6721e-04 | 17.3 | 0.0031 | 6 |
| d15 | 15 | 9 |  | (1,+3) + (1,+4) -> (2,+9) | runaway | 2.47e+08 | 5.5972e-04 | 17.1 | 0.0040 | 6 |
| d17 | 17 | 11 | (2,+1) | (1,+3) + (1,+4) -> (2,+9) | runaway | 3.29e+08 | 4.1212e-04 | 22.6 | 0.0060 | 7 |

Same run from some N on: not reached.

### ladder_2p/l1_n+3_l1_n+4_ladder2

Parents: (1,+3) f 15.974 c/d, gamma -5.538e-09 s^-1; (1,+4) f 17.822 c/d, gamma -7.208e-08 s^-1

Built by: step k: the k daughters shared by both parents with the lowest larger E_th of their two pairs, each with both pairs; plus, for each combination of the parents, its direct triplet of highest R_c.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (1,+3) -> (4,-1) + (5,-58) | +1.163e-02 | +4.294e+00 | 1.483e-08 | none |  |
| pair | (1,+4) -> (3,-15) + (4,-1) | +3.022e-02 | +9.673e-01 | 3.703e-08 | none |  |
| pair | (1,+3) -> (1,-14) + (2,+1) | +2.560e-03 | +4.848e-01 | 3.021e-09 | none |  |
| pair | (1,+4) -> (2,+1) + (3,-15) | -1.627e-02 | +6.737e-01 | 3.939e-08 | none |  |
| pair | (1,+3) -> (2,+0) + (3,-15) | +4.824e-02 | +3.639e-01 | 1.407e-07 | none |  |
| pair | (1,+3) -> (5,-2) + (6,-67) | -2.503e-03 | +1.797e+00 | 3.653e-10 | 3.035e-07 |  |
| pair | (1,+4) -> (5,-2) + (6,-29) | -1.706e-02 | +9.372e-01 | 3.388e-07 | none |  |
| pair | (1,+3) -> (3,-1) + (4,-29) | +2.270e-02 | +1.298e+00 | 4.159e-08 | none |  |
| pair | (1,+4) -> (3,-1) + (4,-15) | +9.330e-02 | +4.754e-01 | 3.746e-07 | none |  |
| pair | (1,+3) -> (4,-15) + (5,-4) | +4.205e-02 | +3.614e-01 | 8.152e-08 | none |  |
| pair | (1,+3) -> (3,+0) + (4,-41) | -4.766e-03 | +2.332e+00 | 9.719e-09 | none |  |
| pair | (1,+4) -> (3,+0) + (4,-18) | +9.113e-02 | +9.054e-01 | 4.604e-07 | none |  |
| pair | (1,+3) -> (4,-18) + (5,-3) | -1.253e-02 | +1.875e-01 | 7.388e-07 | none |  |
| sum | (1,+3) + (1,+4) -> (2,+9) | +1.995e-03 | +8.019e+00 |  |  | (R_c 6.170e+07) |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d4 | 4 | 3 | (4,-1) | (1,+3) + (1,+4) -> (2,+9) | runaway | 1.98e+08 | 6.9560e-04 | 13.7 | 0.0036 | 2 |
| d6 | 6 | 5 | (2,+1) (4,-1) | (1,+3) + (1,+4) -> (2,+9) | runaway | 2.04e+08 | 6.7576e-04 | 14.8 | 0.0023 | 3 |
| d7 | 7 | 6 | (2,+1) (4,-1) (3,-15) | (1,+3) + (1,+4) -> (2,+9) | runaway | 2.04e+08 | 6.7541e-04 | 14.8 | 0.0023 | 3 |
| d10 | 10 | 8 | (2,+1) (4,-1) (5,-2) (3,-15) | (1,+3) + (1,+4) -> (2,+9) | runaway | 2.21e+08 | 6.5096e-04 | 13.9 | 0.0023 | 5 |
| d13 | 13 | 10 | (2,+1) (4,-1) (5,-2) (3,-15) (3,-1) | (1,+3) + (1,+4) -> (2,+9) | runaway | 2.35e+08 | 5.6733e-04 | 17 | 0.0042 | 7 |
| d14 | 14 | 11 | (2,+1) (4,-1) (5,-2) (3,-15) (3,-1) (4,-15) | (1,+3) + (1,+4) -> (2,+9) | runaway | 2.39e+08 | 6.1274e-04 | 15.9 | 0.0029 | 7 |
| d17 | 17 | 13 | (2,+1) (3,+0) (4,-1) (5,-2) (3,-15) (3,-1) (4,-15) | (1,+3) + (1,+4) -> (2,+9) | runaway | 2.44e+08 | 5.6438e-04 | 16.8 | 0.0038 | 9 |
| d18 | 18 | 14 | (2,+1) (3,+0) (4,-1) (5,-2) (3,-15) (3,-1) (4,-15) (4,-18) | (1,+3) + (1,+4) -> (2,+9) | runaway | 2.45e+08 | 5.5142e-04 | 17.3 | 0.0041 | 7 |

Same run from some N on: not reached.

### ladder_2p/l1_n+3_l1_n+5

Parents: (1,+3) f 15.974 c/d, gamma -5.538e-09 s^-1; (1,+5) f 20.671 c/d, gamma -3.430e-07 s^-1

Built by: step k: the k daughters shared by both parents with the lowest larger E_th of their two pairs, each with both pairs.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (1,+3) -> (6,-5) + (7,-29) | -6.572e-03 | +1.945e+00 | 4.755e-10 | none |  |
| pair | (1,+5) -> (5,-8) + (6,-5) | +1.062e-03 | +2.510e+00 | 1.773e-11 | none |  |
| pair | (1,+3) -> (3,+0) + (4,-41) | -4.766e-03 | +2.332e+00 | 9.719e-09 | none |  |
| pair | (1,+5) -> (2,-5) + (3,+0) | -8.364e-03 | +2.064e-01 | 1.661e-08 | none |  |
| pair | (1,+3) -> (5,-4) + (6,-23) | -4.781e-02 | +1.077e+00 | 4.471e-08 | none |  |
| pair | (1,+5) -> (5,-4) + (6,-9) | -3.343e-02 | +3.236e+00 | 2.678e-08 | none |  |
| pair | (1,+3) -> (1,+1) + (2,-8) | -1.742e-02 | +3.586e-02 | 3.318e-07 | none |  |
| pair | (1,+5) -> (1,+1) + (2,-3) | -3.336e-02 | +1.296e+00 | 1.382e-07 | none |  |
| pair | (1,+3) -> (1,+2) + (2,-23) | +1.014e-02 | +8.594e-01 | 1.861e-08 | none |  |
| pair | (1,+5) -> (1,+2) + (2,-5) | +7.114e-02 | +1.951e-01 | 1.174e-06 | none |  |
| pair | (1,+3) -> (3,-9) + (4,-4) | -1.816e-02 | +1.099e-01 | 5.203e-07 | none |  |
| pair | (1,+5) -> (4,-4) + (5,-6) | -1.381e-02 | +1.578e-01 | 1.927e-06 | none |  |
| pair | (1,+3) -> (4,-3) + (5,-17) | +1.255e-01 | +4.480e-01 | 7.304e-07 | none |  |
| pair | (1,+5) -> (3,-4) + (4,-3) | -1.597e-01 | +1.172e+00 | 3.061e-06 | none |  |
| pair | (1,+3) -> (3,-2) + (4,-15) | -1.001e-01 | +1.943e-01 | 3.410e-06 | none |  |
| pair | (1,+5) -> (2,-3) + (3,-2) | -4.559e-02 | +9.813e-01 | 7.993e-07 | none |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d3 | 3 | 2 | (6,-5) |  | runaway | 3.84e+07 | 7.6550e-04 | 12.8 | 0.0027 | 2 |
| d6 | 6 | 4 | (3,+0) (6,-5) |  | runaway | 3.84e+07 | 7.6545e-04 | 12.8 | 0.0027 | 2 |
| d9 | 9 | 6 | (3,+0) (6,-5) (5,-4) |  | runaway | 4.27e+07 | 6.8679e-04 | 14.5 | 0.0027 | 4 |
| d12 | 12 | 8 | (3,+0) (6,-5) (5,-4) (1,+1) |  | runaway | 4.27e+07 | 6.9745e-04 | 14.1 | 0.0025 | 4 |
| d14 | 14 | 10 | (1,+2) (3,+0) (6,-5) (5,-4) (1,+1) |  | runaway | 4.28e+07 | 7.0314e-04 | 13.9 | 0.0026 | 4 |
| d17 | 17 | 12 | (1,+2) (3,+0) (6,-5) (5,-4) (1,+1) (4,-4) |  | runaway | 4.29e+07 | 6.7701e-04 | 14.8 | 0.0029 | 4 |
| d20 | 20 | 14 | (1,+2) (3,+0) (6,-5) (5,-4) (1,+1) (4,-4) (4,-3) |  | runaway | 4.81e+07 | 6.1993e-04 | 15.4 | 0.0037 | 6 |
| d22 | 22 | 16 | (1,+2) (3,+0) (6,-5) (5,-4) (1,+1) (4,-4) (4,-3) (3,-2) |  | runaway | 4.48e+07 | 6.5545e-04 | 14.8 | 0.0029 | 6 |

Same run from some N on: not reached.

### ladder_2p/l1_n+3_l1_n+5_ladder2

Parents: (1,+3) f 15.974 c/d, gamma -5.538e-09 s^-1; (1,+5) f 20.671 c/d, gamma -3.430e-07 s^-1

Built by: step k: the k daughters shared by both parents with the lowest larger E_th of their two pairs, each with both pairs; plus, for each combination of the parents, its direct triplet of highest R_c.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (1,+3) -> (6,-5) + (7,-29) | -6.572e-03 | +1.945e+00 | 4.755e-10 | none |  |
| pair | (1,+5) -> (5,-8) + (6,-5) | +1.062e-03 | +2.510e+00 | 1.773e-11 | none |  |
| pair | (1,+3) -> (3,+0) + (4,-41) | -4.766e-03 | +2.332e+00 | 9.719e-09 | none |  |
| pair | (1,+5) -> (2,-5) + (3,+0) | -8.364e-03 | +2.064e-01 | 1.661e-08 | none |  |
| pair | (1,+3) -> (5,-4) + (6,-23) | -4.781e-02 | +1.077e+00 | 4.471e-08 | none |  |
| pair | (1,+5) -> (5,-4) + (6,-9) | -3.343e-02 | +3.236e+00 | 2.678e-08 | none |  |
| pair | (1,+3) -> (1,+1) + (2,-8) | -1.742e-02 | +3.586e-02 | 3.318e-07 | none |  |
| pair | (1,+5) -> (1,+1) + (2,-3) | -3.336e-02 | +1.296e+00 | 1.382e-07 | none |  |
| pair | (1,+3) -> (1,+2) + (2,-23) | +1.014e-02 | +8.594e-01 | 1.861e-08 | none |  |
| pair | (1,+5) -> (1,+2) + (2,-5) | +7.114e-02 | +1.951e-01 | 1.174e-06 | none |  |
| pair | (1,+3) -> (3,-9) + (4,-4) | -1.816e-02 | +1.099e-01 | 5.203e-07 | none |  |
| pair | (1,+5) -> (4,-4) + (5,-6) | -1.381e-02 | +1.578e-01 | 1.927e-06 | none |  |
| pair | (1,+3) -> (4,-3) + (5,-17) | +1.255e-01 | +4.480e-01 | 7.304e-07 | none |  |
| pair | (1,+5) -> (3,-4) + (4,-3) | -1.597e-01 | +1.172e+00 | 3.061e-06 | none |  |
| pair | (1,+3) -> (3,-2) + (4,-15) | -1.001e-01 | +1.943e-01 | 3.410e-06 | none |  |
| pair | (1,+5) -> (2,-3) + (3,-2) | -4.559e-02 | +9.813e-01 | 7.993e-07 | none |  |
| sum | (1,+3) + (1,+5) -> (2,+10) | +2.383e-01 | +1.364e+01 |  |  | (R_c 1.389e+07) |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d4 | 4 | 3 | (6,-5) | (1,+3) + (1,+5) -> (2,+10) | runaway | 3.84e+07 | 7.7085e-04 | 13 | 0.0026 | 2 |
| d7 | 7 | 5 | (3,+0) (6,-5) | (1,+3) + (1,+5) -> (2,+10) | runaway | 3.84e+07 | 7.7121e-04 | 13 | 0.0026 | 2 |
| d10 | 10 | 7 | (3,+0) (6,-5) (5,-4) | (1,+3) + (1,+5) -> (2,+10) | runaway | 4.28e+07 | 6.9376e-04 | 14.3 | 0.0030 | 4 |
| d13 | 13 | 9 | (3,+0) (6,-5) (5,-4) (1,+1) | (1,+3) + (1,+5) -> (2,+10) | runaway | 4.29e+07 | 6.9846e-04 | 14.2 | 0.0028 | 4 |
| d15 | 15 | 11 | (1,+2) (3,+0) (6,-5) (5,-4) (1,+1) | (1,+3) + (1,+5) -> (2,+10) | runaway | 4.25e+07 | 6.8174e-04 | 14.4 | 0.0030 | 4 |
| d18 | 18 | 13 | (1,+2) (3,+0) (6,-5) (5,-4) (1,+1) (4,-4) | (1,+3) + (1,+5) -> (2,+10) | runaway | 4.27e+07 | 6.9699e-04 | 14.2 | 0.0024 | 4 |
| d21 | 21 | 15 | (1,+2) (3,+0) (6,-5) (5,-4) (1,+1) (4,-4) (4,-3) | (1,+3) + (1,+5) -> (2,+10) | runaway | 4.58e+07 | 6.4537e-04 | 15 | 0.0035 | 6 |
| d23 | 23 | 17 | (1,+2) (3,+0) (6,-5) (5,-4) (1,+1) (4,-4) (4,-3) (3,-2) | (1,+3) + (1,+5) -> (2,+10) | runaway | 4.87e+07 | 6.1326e-04 | 15.6 | 0.0044 | 6 |

Same run from some N on: not reached.

