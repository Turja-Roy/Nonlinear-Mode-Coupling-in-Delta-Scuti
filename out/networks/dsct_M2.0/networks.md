# dsct_M2.0: ladders and runs

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

### ladder_1p/l1_n-1

Parents: (1,-1) f 10.424 c/d, gamma -2.276e-10 s^-1

Built by: pairs ranked by E_th, lowest first; a pair that would bring N past the target skipped.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (1,-1) -> (4,-38) + (5,-9) | +7.011e-04 | +3.417e-01 | 1.333e-10 | none |  |
| pair | (1,-1) -> (4,-8) + (5,-37) | -1.285e-03 | +1.712e-01 | 2.459e-09 | none |  |
| pair | (1,-1) -> (5,-33) + (6,-12) | -1.096e-03 | +2.071e-01 | 9.314e-09 | none |  |
| pair | (1,-1) -> (6,-14) + (7,-31) | -7.146e-04 | +1.617e-01 | 1.054e-08 | none |  |
| pair | (1,-1) -> (5,-9) + (6,-56) | -7.120e-03 | +8.545e-02 | 1.750e-08 | 2.071e-08 |  |
| pair | (1,-1) -> (5,-9) + (6,-55) | +1.924e-02 | +8.889e-02 | 1.287e-07 | 1.551e-07 |  |
| pair | (1,-1) -> (6,-11) + (7,-54) | +7.867e-03 | +4.988e-02 | 1.871e-07 | 2.195e-07 |  |
| pair | (1,-1) -> (5,-9) + (6,-57) | -4.049e-02 | +9.321e-02 | 3.273e-07 | 3.662e-07 |  |
| pair | (1,-1) -> (5,-39) + (6,-11) | +1.049e-02 | +1.493e-01 | 4.544e-07 | 1.624e-04 |  |
| pair | (1,-1) -> (5,-9) + (6,-54) | +5.458e-02 | +1.072e-01 | 5.942e-07 | 6.958e-07 |  |
| pair | (1,-1) -> (4,-39) + (5,-9) | -5.337e-02 | +3.177e-01 | 8.307e-07 | 1.373e-04 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d2 | 2 | 1 |  |  | runaway | 1.32e+12 | 4.5293e-05 | 4.03 | 0.3434 | 2 |
| d4 | 4 | 2 |  |  | end | 1e+13 | 5.2519e-10 | 0.118 | 1.0006 | 2 |
| d6 | 6 | 3 |  |  | end | 1e+13 | 5.2556e-10 | 0.014 | 1.0000 | 2 |
| d8 | 8 | 4 |  |  | end | 1e+13 | 5.1108e-10 | 4.1 | 1.0247 | 4 |
| d10 | 10 | 6 |  |  | end | 1e+13 | 5.1902e-10 | 3.07 | 0.9932 | 4 |
| d12 | 12 | 7 |  |  | end | 1e+13 | 5.1902e-10 | 3.07 | 0.9932 | 4 |
| d14 | 14 | 9 |  |  | end | 1e+13 | 5.1035e-10 | 3.1 | 1.0008 | 4 |
| d16 | 16 | 11 |  |  | end | 1e+13 | 4.9711e-10 | 3.17 | 1.0144 | 4 |

Same run from some N on: not reached.

## ladder_2p

### ladder_2p/l2_n-2_l1_n-1

Parents: (1,-1) f 10.424 c/d, gamma -2.276e-10 s^-1; (2,-2) f 10.602 c/d, gamma -1.619e-10 s^-1

Built by: step k: the k daughters shared by both parents with the lowest larger E_th of their two pairs, each with both pairs.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (1,-1) -> (5,-33) + (6,-12) | -1.096e-03 | +2.071e-01 | 9.314e-09 | none |  |
| pair | (2,-2) -> (4,-25) + (6,-12) | +9.322e-04 | +2.520e-01 | 1.448e-08 | none |  |
| pair | (1,-1) -> (4,-38) + (5,-9) | +7.011e-04 | +3.417e-01 | 1.333e-10 | none |  |
| pair | (2,-2) -> (5,-9) + (7,-59) | -1.152e-02 | +4.778e-02 | 9.806e-08 | 1.068e-07 |  |
| pair | (1,-1) -> (6,-11) + (7,-54) | +7.867e-03 | +4.988e-02 | 1.871e-07 | 2.195e-07 |  |
| pair | (2,-2) -> (6,-11) + (8,-57) | +1.438e-03 | +2.299e-02 | 1.965e-08 | 2.128e-08 |  |
| pair | (1,-1) -> (4,-8) + (5,-37) | -1.285e-03 | +1.712e-01 | 2.459e-09 | none |  |
| pair | (2,-2) -> (4,-8) + (6,-41) | -1.429e-02 | +9.775e-02 | 4.188e-07 | 1.909e-06 |  |
| pair | (1,-1) -> (5,-9) + (6,-57) | -4.049e-02 | +9.321e-02 | 3.273e-07 | 3.662e-07 |  |
| pair | (2,-2) -> (6,-57) + (6,-10) | +2.056e-02 | +5.107e-02 | 6.966e-07 | 7.540e-07 |  |
| pair | (1,-1) -> (4,-39) + (5,-9) | -5.337e-02 | +3.177e-01 | 8.307e-07 | 1.373e-04 |  |
| pair | (2,-2) -> (4,-39) + (6,-10) | +7.678e-03 | +2.423e-01 | 7.283e-08 | 6.093e-07 |  |
| pair | (1,-1) -> (5,-9) + (6,-58) | -7.655e-02 | +8.813e-02 | 1.185e-06 | 1.309e-06 |  |
| pair | (2,-2) -> (6,-58) + (6,-10) | -1.550e-02 | +4.915e-02 | 3.870e-07 | 4.151e-07 |  |
| pair | (1,-1) -> (6,-10) + (7,-72) | +4.425e-02 | +3.459e-02 | 2.094e-06 | 2.158e-06 |  |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d3 | 3 | 2 | (6,-12) |  | runaway | 8.92e+10 | 5.5071e-04 | 14 | 0.0649 | 3 |
| d6 | 6 | 4 | (5,-9) (6,-12) |  | runaway | 1.91e+11 | 3.5130e-04 | 8.99 | 0.0804 | 3 |
| d9 | 9 | 6 | (5,-9) (6,-12) (6,-11) |  | runaway | 1.91e+11 | 3.5115e-04 | 9 | 0.0804 | 3 |
| d12 | 12 | 8 | (5,-9) (4,-8) (6,-12) (6,-11) |  | runaway | 1.92e+11 | 3.5083e-04 | 8.77 | 0.0805 | 3 |
| d14 | 14 | 10 | (5,-9) (4,-8) (6,-12) (6,-11) (6,-57) |  | end | 1e+13 | 2.5230e-07 | 11 | 1.0116 | 4 |
| d15 | 15 | 12 | (5,-9) (4,-8) (6,-12) (6,-11) (4,-39) (6,-57) |  | end | 1e+13 | 1.2671e-07 | 0.125 | 1.0044 | 1 |
| d16 | 16 | 14 | (5,-9) (4,-8) (6,-12) (6,-11) (4,-39) (6,-57) (6,-58) |  | end | 1e+13 | 6.8511e-08 | 0.371 | 0.9996 | 2 |
| d17 | 17 | 15 | (5,-9) (4,-8) (6,-12) (6,-11) (4,-39) (6,-10) (6,-57) (6,-58) |  | end | 1e+13 | 6.7860e-08 | 0.349 | 1.0014 | 3 |

Same run from some N on: not reached.

### ladder_2p/l2_n-2_l2_n+5

Parents: (2,-2) f 10.602 c/d, gamma -1.619e-10 s^-1; (2,+5) f 28.426 c/d, gamma -4.771e-07 s^-1

Built by: each parent's pairs ranked by E_th, the parents in turn; plus, for each combination of the parents, its direct triplet of highest R_c.

Selecting triplets of the largest network, in the order of its header:

| kind | triplet | Delta [c/d] | kappa | E_th [E_star] | E_eq [E_star] | mu |
|---|---|---|---|---|---|---|
| pair | (2,-2) -> (1,-3) + (3,-11) | -2.977e-04 | +1.407e-01 | 2.759e-09 | none |  |
| pair | (2,-2) -> (4,-25) + (6,-12) | +9.322e-04 | +2.520e-01 | 1.448e-08 | none |  |
| pair | (2,-2) -> (6,-11) + (8,-57) | +1.438e-03 | +2.299e-02 | 1.965e-08 | 2.128e-08 |  |
| pair | (2,-2) -> (4,-39) + (6,-10) | +7.678e-03 | +2.423e-01 | 7.283e-08 | 6.093e-07 |  |
| pair | (2,-2) -> (6,-43) + (6,-11) | +2.265e-03 | +5.481e-02 | 7.549e-08 | 1.768e-07 |  |
| pair | (2,-2) -> (5,-9) + (7,-59) | -1.152e-02 | +4.778e-02 | 9.806e-08 | 1.068e-07 |  |
| pair | (2,-2) -> (6,-10) + (8,-76) | -7.988e-03 | +2.196e-02 | 1.008e-07 | 1.022e-07 |  |
| pair | (2,-2) -> (5,-9) + (7,-58) | +2.300e-02 | +6.139e-02 | 1.926e-07 | 2.067e-07 |  |
| pair | (2,-2) -> (6,-10) + (8,-75) | +1.238e-02 | +1.997e-02 | 3.677e-07 | 3.741e-07 |  |
| pair | (2,-2) -> (6,-58) + (6,-10) | -1.550e-02 | +4.915e-02 | 3.870e-07 | 4.151e-07 |  |
| pair | (2,-2) -> (5,-43) + (5,-9) | -2.246e-02 | +1.170e-01 | 5.225e-07 | 1.241e-06 |  |
| sum | (2,-2) + (2,+5) -> (2,+9) | -9.788e-03 | +2.971e+00 |  |  | (R_c 2.415e+07) |
| harmonic | (2,+5) + (2,+5) -> (2,+14) | +3.267e-01 | +4.369e+01 |  |  | (R_c 3.365e+07) |

| run | N | triplets | shared daughters | direct triplets | stop | t_stop [s] | <E_par> [E_star] | swing | dissipation/driving | n_active |
|---|---|---|---|---|---|---|---|---|---|---|
| d4 | 4 | 3 |  | (2,-2) + (2,+5) -> (2,+9); (2,+5) + (2,+5) -> (2,+14) | end | 1e+13 | 7.3937e-10 | 1.22e-10 | 1.0057 | 1 |
| d6 | 6 | 4 |  | (2,-2) + (2,+5) -> (2,+9); (2,+5) + (2,+5) -> (2,+14) | end | 1e+13 | 7.3937e-10 | 1.19e-10 | 1.0057 | 1 |
| d8 | 8 | 5 |  | (2,-2) + (2,+5) -> (2,+9); (2,+5) + (2,+5) -> (2,+14) | end | 1e+13 | 7.3937e-10 | 1.18e-10 | 1.0057 | 1 |
| d10 | 10 | 6 |  | (2,-2) + (2,+5) -> (2,+9); (2,+5) + (2,+5) -> (2,+14) | end | 1e+13 | 7.3937e-10 | 1.17e-10 | 1.0057 | 1 |
| d12 | 12 | 8 |  | (2,-2) + (2,+5) -> (2,+9); (2,+5) + (2,+5) -> (2,+14) | end | 1e+13 | 7.3937e-10 | 1.2e-10 | 1.0057 | 1 |
| d14 | 14 | 9 |  | (2,-2) + (2,+5) -> (2,+9); (2,+5) + (2,+5) -> (2,+14) | end | 1e+13 | 7.3937e-10 | 1.21e-10 | 1.0057 | 1 |
| d16 | 16 | 11 |  | (2,-2) + (2,+5) -> (2,+9); (2,+5) + (2,+5) -> (2,+14) | end | 1e+13 | 7.3937e-10 | 1.2e-10 | 1.0057 | 1 |
| d18 | 18 | 13 |  | (2,-2) + (2,+5) -> (2,+9); (2,+5) + (2,+5) -> (2,+14) | end | 1e+13 | 7.3937e-10 | 2.21e-09 | 1.0057 | 1 |

Same run (stop, n_active, swing within 0.01, <E_par> within 1e-3) from N = 4 on.

