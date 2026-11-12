# dsct_M1.7_T7750: ladders and runs

Written by scripts/network_index.py.

Notation. (l,n): degree and radial order of a mode; f: frequency; gamma: linear rate, < 0 driven (parent),
> 0 damped (daughter). N: number of daughters. Pair p -> b + c: parametric triplet; harmonic x + x -> c and
sum x + y -> c: direct triplets. Delta: detuning; kappa: coupling coefficient. E_th: parent energy above
which a pair grows (Weinberg+12 eq. 50); E_eq: parent energy at a triplet's fixed point (MW25 eq. 6), for a
sum triplet E_x + E_y. Energies in E_star = GM^2/R. Shared daughter: coupled to both parents. Run columns,
second half of the run: <E_par>, mean summed parent energy; swing, (max - min)/mean of it;
dissipation/driving, sum 2 gamma_d <E_d> over sum 2|gamma_p| <E_p>; n_active, daughters above 1e-3 of the
largest. stop: end (t = 1e13 s), runaway (summed energy 1e-2), wall (time limit).

## ladder_2p

### ladder_2p/l0_n+1_l2_n+6

Parents: (0,+1) f 18.798 c/d, gamma -8.280e-11 s^-1; (2,+6) f 52.612 c/d, gamma -2.985e-08 s^-1

Built by: each parent's pairs ranked by E_th, the parents in turn; plus, for each combination of the parents, its direct triplet of highest R_c.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (0,+1) -> (6,-26) + (6,-6) | +3.026e-02 | +8.174e-02 | 3.618e-07 | 4.744e-07 |  |
| pair | (0,+1) -> (4,-12) + (4,-5) | +5.758e-02 | +3.330e-01 | 3.134e-06 | none |  |
| pair | (0,+1) -> (6,-18) + (6,-7) | +4.893e-02 | +2.054e-01 | 3.939e-06 | none |  |
| pair | (0,+1) -> (6,-27) + (6,-6) | -1.081e-01 | +6.373e-02 | 5.814e-06 | 7.078e-06 |  |
| pair | (0,+1) -> (4,-9) + (4,-6) | -1.691e-02 | +1.370e-01 | 7.660e-06 | none |  |
| pair | (0,+1) -> (5,-15) + (5,-6) | +7.410e-02 | +2.702e-01 | 8.406e-06 | none |  |
| pair | (0,+1) -> (6,-25) + (6,-6) | +1.831e-01 | +8.612e-02 | 1.621e-05 | 2.403e-05 |  |
| pair | (0,+1) -> (6,-28) + (6,-6) | -2.304e-01 | +6.977e-02 | 1.813e-05 | 2.118e-05 |  |
| pair | (0,+1) -> (6,-30) + (6,-6) | -4.370e-01 | +6.182e-02 | 5.419e-05 | 5.955e-05 |  |
| pair | (0,+1) -> (6,-24) + (6,-6) | +3.499e-01 | +9.928e-02 | 6.147e-05 | 1.111e-04 |  |
| pair | (0,+1) -> (6,-19) + (6,-7) | -2.328e-01 | +1.845e-01 | 8.198e-05 | none |  |
| sum | (0,+1) + (2,+6) -> (2,+9) | -1.996e-01 | +3.396e+00 |  |  | (R_c 5.566e+06) |
| harmonic | (2,+6) + (2,+6) -> (0,+15) | -5.470e-01 | +1.041e+02 |  |  | (R_c 3.029e+08) |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d4 | 4 | 3 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 4.3336e-07 | 44.6 | 1.0294 | 1 |
| d6 | 6 | 4 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | runaway | 1.98e+11 | 6.2448e-04 | 6.12 | 0.0253 | 2 |
| d8 | 8 | 5 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | runaway | 1.98e+11 | 6.2560e-04 | 6.04 | 0.0253 | 2 |
| d10 | 10 | 7 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 4.0230e-07 | 124 | 1.4459 | 2 |
| d12 | 12 | 8 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 4.1114e-07 | 17.6 | 0.9947 | 2 |
| d14 | 14 | 9 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | end | 1e+13 | 4.6136e-07 | 15.1 | 0.9339 | 3 |
| d16 | 16 | 11 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | wall | 5.2e+12 | 4.5049e-07 | 27.5 | 0.8064 | 3 |
| d18 | 18 | 13 |  | (0,+1) + (2,+6) -> (2,+9); (2,+6) + (2,+6) -> (0,+15) | wall | 5.18e+12 | 4.6267e-07 | 27.6 | 0.7669 | 3 |

Same run from some N on: not reached.

### ladder_2p/l2_n-3_l0_n+7

Parents: (0,+7) f 53.136 c/d, gamma -1.697e-08 s^-1; (2,-3) f 11.802 c/d, gamma -3.458e-14 s^-1

Built by: each parent's pairs ranked by E_th, the parents in turn; plus, for each combination of the parents, its direct triplet of highest R_c.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (2,-3) -> (6,-27) + (8,-17) | +5.995e-04 | +1.564e-01 | 3.752e-09 | 3.752e-09 |  |
| pair | (2,-3) -> (5,-27) + (5,-10) | -3.930e-04 | -2.191e-02 | 7.260e-09 | 7.261e-09 |  |
| pair | (2,-3) -> (11,-43) + (13,-30) | +9.303e-04 | -1.596e-01 | 1.679e-08 | 1.679e-08 |  |
| pair | (2,-3) -> (2,-19) + (4,-7) | +1.453e-03 | +9.292e-02 | 2.952e-08 | 2.955e-08 |  |
| pair | (2,-3) -> (11,-54) + (13,-26) | +1.392e-03 | +6.656e-02 | 4.293e-08 | 4.293e-08 |  |
| pair | (2,-3) -> (13,-70) + (15,-29) | -7.609e-04 | +1.510e-01 | 4.862e-08 | 4.862e-08 |  |
| pair | (2,-3) -> (12,-56) + (14,-29) | +1.296e-03 | -6.429e-02 | 9.344e-08 | 9.344e-08 |  |
| pair | (2,-3) -> (14,-37) + (16,-52) | +1.438e-03 | +1.336e-01 | 1.130e-07 | 1.130e-07 |  |
| sum | (0,+7) + (2,-3) -> (2,+8) | -1.107e-02 | +2.218e-01 |  |  | (R_c 1.583e+06) |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d3 | 3 | 2 |  | (0,+7) + (2,-3) -> (2,+8) | runaway | 6.79e+08 | 8.7155e-04 | 11.5 | 0.0000 | 3 |
| d5 | 5 | 3 |  | (0,+7) + (2,-3) -> (2,+8) | runaway | 6.79e+08 | 8.7160e-04 | 11.5 | 0.0000 | 5 |
| d7 | 7 | 4 |  | (0,+7) + (2,-3) -> (2,+8) | runaway | 6.79e+08 | 8.7156e-04 | 11.5 | 0.0000 | 6 |
| d9 | 9 | 5 |  | (0,+7) + (2,-3) -> (2,+8) | runaway | 6.79e+08 | 8.7157e-04 | 11.5 | 0.0000 | 8 |
| d11 | 11 | 6 |  | (0,+7) + (2,-3) -> (2,+8) | runaway | 6.79e+08 | 8.7162e-04 | 11.5 | 0.0000 | 9 |
| d13 | 13 | 7 |  | (0,+7) + (2,-3) -> (2,+8) | runaway | 6.79e+08 | 8.7164e-04 | 11.5 | 0.0000 | 10 |
| d15 | 15 | 8 |  | (0,+7) + (2,-3) -> (2,+8) | runaway | 6.79e+08 | 8.7162e-04 | 11.5 | 0.0000 | 11 |
| d17 | 17 | 9 |  | (0,+7) + (2,-3) -> (2,+8) | runaway | 6.79e+08 | 8.7156e-04 | 11.5 | 0.0000 | 12 |

Same run from some N on: not reached.

