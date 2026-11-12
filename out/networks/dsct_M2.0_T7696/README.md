# Mode networks of dsct_M2.0_T7696

> Folders renamed 2026-11-11 (after this file was written); `networks.md` beside it lists every run.
> `ladder_1p/cut0.15_d<N>` → `ladder_1p/l2_n-2/d<N>`;
> `ladder_2p_shared/cut0.15_shared<k>_d<N>` → `ladder_2p/l2_n-3_l2_n-2/d<N>`;
> `ladder_2p_direct/cut0.15_d<N>_dir1` → `ladder_2p/l2_n-2_l1_n+7/d<N>`;
> `_noclosure`, `ladder_1p_direct`, `ladder_2p_shared_direct`: deleted, the same networks as the sets above.

2 M☉, log g = 3.93, T_eff = 7695 K (MW23's representative model). Written 2026-11-11 from the
network files and run summaries in this folder; nothing was rerun for it.

## Notation

- (l, n): spherical degree l and radial order n of a GYRE mode; n < 0 is a g mode, n > 0 a p mode.
  Every mode here has m = 0.
- f [c/d]: mode frequency in cycles per day; ω = 2πf.
- γ [s⁻¹]: linear rate, γ = γ_rad + γ_turb (radiative + turbulent). γ < 0: the mode is driven
  (self-excited); γ > 0: damped.
- Parent: a driven mode. Daughter: a damped mode.
- E: mode energy in units of E★ = GM²/R.
- κ: three-mode coupling coefficient (dimensionless).
- Δ [c/d]: detuning, how far a triplet is from exact resonance.
  Parametric triplet a → b + c (parent a feeds daughters b and c): Δ = f_b + f_c − f_a.
  Direct triplet x + y → c (parents x and y feed a damped mode c near f_x + f_y): Δ = f_c − f_x − f_y.
  x = y is allowed (x + x → c, c near the harmonic 2f_x).
- E_th: parametric threshold. The pair b + c grows once the parent's energy exceeds E_th (Weinberg+12 eq. 50).
- E_eq (parametric): parent energy at the stable fixed point of a → b + c (MW25 eq. 6). "none" means
  γ_a + γ_b + γ_c ≤ 0, so no stable fixed point exists.
- E_eq (direct): E_x + E_y at the fixed point of x + y → c. For x = y it is E_x, the energy at which the
  triplet holds x.
- μ = |ω_c κ| / √(Δ² + γ_c²): direct-coupling strength, |A_c| = μ |A_x||A_y| (MW23 eq. 9).
- R_c = (sμ)², with s = 2 for x ≠ y and 1 for x = y: the direct-triplet ranking used before 2026-11-10.
- Run outcome, taken over the second half of each run:
  - ⟨E_parents⟩: mean total parent energy.
  - swing: (max − min)/mean of the total parent energy; 0 means steady.
  - dissipation/driving: mean damping power over mean driving power, with P = Σ 2|γ_p|E_p and D = Σ 2γ_d E_d
    (bottom panel of every plot). 1 means driving and damping balance.
  - active daughters: daughters whose mean energy is above 10⁻³ of the largest daughter's.
  - runaway: the run stopped because the total energy passed 10⁻² E★.

## What is in this folder

| folder | parents | one step adds | networks | direct triplets |
|---|---|---|---|---|
| `ladder_1p` | (2,−2) | the next parametric pair of (2,−2) | d2 … d16, each also without closure (`_noclosure`) | none |
| `ladder_1p_direct` | (2,−2) | as `ladder_1p`, plus the best x + x → c | d2_dir1 … d16_dir1 | **none exist** (see Updates 2) |
| `ladder_2p_shared` | (2,−2), (2,−3) | one shared daughter with its two pairs | shared1_d3 … shared8_d21 | none |
| `ladder_2p_shared_direct` | (2,−2), (2,−3) | as `ladder_2p_shared`, plus the best direct triplets | shared1_d3_dir1 … shared8_d21_dir1 | **none exist** (see Updates 2) |
| `ladder_2p_direct` | (1,+7), (2,−2) | the next parametric pair | d4_dir1 … d18_dir1 | 2: the mixed networks (MW25 sec. 5) |

Folder names:
- `cut0.15`: only triplets with |Δ| < 0.15 √(GM/R³) = 0.453 c/d are candidates.
- `d<N>`: N damped modes in the network.
- `shared<k>`: k shared daughters.
- `_dir1`: the best direct triplet for each pair of parents was added.

Each run folder holds:
- `network.data`: modes and triplets, with the selecting pairs in the header.
- `energies.csv`: E(t) of every mode.
- `summary.csv`: the outcome numbers above.
- `log.txt`
- `network_run.png` / `.pdf`

## How the networks are built

1. **Mode pool.** Every GYRE mode of the model with l ≤ 25, using γ = γ_rad + γ_turb.
2. **Parents.**
   - `ladder_1p` and `ladder_2p_shared`: (2,−2) and (2,−3), the driven l ≤ 2 modes with the lowest E_th
     (the October default, kept fixed since).
   - `ladder_2p_direct`: picked by hand from the direct-triplet scan (next section).
3. **Parametric pairs.**
   - Candidates are parent → b + c with b and c both damped and |Δ| below the cut.
   - Each parent's pairs are ranked by E_th, lowest first.
   - A step adds the next pair. A pair that would take the daughter count past N is skipped.
   - With two parents, the parents take turns.
4. **Shared daughters** (`ladder_2p_shared`).
   - A damped mode d qualifies when both parents have a pair containing it, p₁ → d + x and p₂ → d + y.
   - Qualifying modes are ranked by the larger of the two E_th, since d links the parents only once both
     pairs grow.
   - A step adds d with both of its pairs.
5. **Direct triplets** (`--direct 1`).
   - For each pair of parents, including a parent with itself, the best x + y → c is added, with c damped
     and |Δ| below the cut.
   - Ranked by E_eq today; these files were ranked by R_c. Both give the same picks here (Updates 4).
6. **Closure.** Every further triplet among the chosen modes with |Δ| < 0.005 √(GM/R³) = 0.015 c/d is
   added. It adds none in any network here.
7. **Run settings.**
   - dA_dt, rtol 10⁻⁸, t_end = 10¹³ s.
   - Each parent starts at E = 10⁻¹², each daughter at 10⁻¹⁵, with random phases (seed 1).
   - A run stops early if the total energy passes 10⁻² E★.

## The mixed networks: `ladder_2p_direct`

### Why these two parents

The scan of every direct triplet x + y → c with x and y driven (l ≤ 2) and c damped found 111 of them.
In none of the 111 do both parents have a parametric pair with an E_eq.

- The slow modes that do have such pairs sit at 8–19 c/d. Any two of them add up to a frequency that is
  still in the driven band, so the mode c there is not damped.
- A damped c therefore needs a fast, high-order p-mode partner. Those partners have no damped pairs,
  because half their frequency lies in the driven band.

So the pick is a slow parent with pairs, (2,−2), plus a fast partner, (1,+7), with:
- a sum triplet linking the two: (1,+7) + (2,−2) → (1,+10);
- a harmonic triplet (1,+7) + (1,+7) → (0,+16) to hold (1,+7), since it has no pairs of its own. That
  hold is stable when ζ = |γ_x|/γ_c < 1/2; here ζ = 0.032.

| parent | f [c/d] | γ [s⁻¹] | damped pairs within the cut (with E_eq) |
|---|---|---|---|
| (2,−2) | 10.944 | −6.697e-11 | 505 (247) |
| (1,+7) | 32.096 | −7.808e-07 | 0 |

### Direct triplets (in every network of the set)

| triplet | f_c [c/d] | γ_c [s⁻¹] | Δ [c/d] | κ | E_eq | μ | R_c | other candidates |
|---|---|---|---|---|---|---|---|---|
| (1,+7) + (1,+7) → (0,+16) | 64.200 | 2.454e-05 | +7.40e-3 | 114.7 | 6.687e-11 | 2.181e4 | 4.758e8 | (2,+15): E_eq 1.305e-10, μ 1.612e4 |
| (1,+7) + (2,−2) → (1,+10) | 43.023 | 4.443e-06 | −1.785e-2 | 4.909 | 5.549e-09 | 3.318e3 | 4.404e7 | none |

(2,−2) + (2,−2) has no damped c within the cut.

### Parametric pairs, all of parent (2,−2), in the order they enter

| # | pair | f_b, f_c [c/d] | γ_b, γ_c [s⁻¹] | Δ [c/d] | κ | E_th | E_eq | first in |
|---|---|---|---|---|---|---|---|---|
| 1 | (3,−31) + (5,−8) | 1.769, 9.175 | 8.10e-11, 7.37e-13 | −4.04e-5 | 0.215 | 4.87e-12 | 1.50e-10 | d4 |
| 2 | (5,−9) + (7,−50) | 8.548, 2.397 | 1.11e-12, 2.55e-09 | +8.22e-4 | 0.053 | 1.30e-09 | 1.37e-09 | d6 |
| 3 | (4,−29) + (4,−7) | 2.423, 8.519 | 8.66e-11, 6.83e-13 | −1.66e-3 | 0.175 | 8.45e-09 | 1.56e-07 | d8 |
| 4 | (4,−24) + (6,−11) | 2.906, 8.037 | 2.81e-11, 2.15e-12 | −1.13e-3 | 0.287 | 1.09e-08 | none | d10 |
| 5 | (5,−21) + (5,−11) | 3.952, 6.992 | 1.41e-11, 2.05e-12 | −8.88e-4 | 0.105 | 7.17e-08 | none | d12 |
| 6 | (6,−29) + (6,−12) | 3.470, 7.477 | 7.94e-11, 2.30e-12 | +2.33e-3 | 0.090 | 1.76e-07 | 5.42e-06 | d14 |
| 7 | (6,−49) + (6,−10) | 2.129, 8.823 | 1.95e-09, 1.62e-12 | +7.68e-3 | 0.050 | 2.56e-07 | 2.74e-07 | d16 |
| 8 | (6,−35) + (6,−11)* | 2.911, 8.037 | 1.56e-10, 2.15e-12 | +3.39e-3 | 0.068 | 3.60e-07 | 1.09e-06 | d18 |
| 9 | (6,−10)* + (8,−64) | 8.823, 2.135 | 1.62e-12, 1.13e-08 | +1.33e-2 | 0.021 | 7.87e-07 | 7.97e-07 | d18 |

\* daughter already in the network from an earlier pair. So d18 has 16 parametric daughters plus (0,+16)
and (1,+10), 20 modes and 11 triplets in all.

### What the runs do

All eight networks (d4 … d18) give the same run:

| ⟨E_parents⟩ | swing | dissipation/driving | active daughters | stop |
|---|---|---|---|---|
| 6.687e-11 | 0.000 | 1.0001 | 1, (0,+16) | end, t = 10¹³ s |

- **(1,+7)** grows and settles at 6.687e-11, the E_eq of its harmonic triplet. Its driving,
  2|γ|E = 1.04e-16 E★/s, is all dissipated by (0,+16), which takes 100% of ⟨D⟩.
- **(2,−2)** is drained from the start. It falls below 10⁻²² E★ by t = 2.3 × 10¹⁰ s and below 10⁻³⁴ by
  7.7 × 10¹¹ s.
- **Parametric pairs.** Every pair hangs on (2,−2), so none ever grows. Adding pairs (d4 → d18) changes
  nothing.
- So these networks have only one parent alive at the end, and are not mixed states (MW25 sec. 5 needs
  both parents alive).

### Why (2,−2) dies

Through the sum triplet, (1,+7) drains (2,−2) at a rate proportional to E₍₁,₊₇₎. Treating (1,+10) as
following its forcing:

dE_s/dt = 2E_s ( |γ_s| − γ_c R_c (ω_s/ω_c) E_f ),

where s is the slow parent (2,−2), f the fast parent (1,+7), and c is (1,+10).

The slow parent therefore decays once E_f > K_s = |γ_s| ω_c / (γ_c R_c ω_s) = 1.35e-12. The fast parent is
held at 6.69e-11, about 50 K_s. The same estimate explained the other models' mixed runs (Plans,
network-growth.md); a stable mixed state needs the fast parent held below K_s. That has not been tested.

## `ladder_1p` (and `ladder_1p_direct`)

- The parent is (2,−2), with the pairs of the table above in the same order.
- `ladder_1p` d2 holds pair 1, d4 adds pair 2, and so on up to d14 (pair 7). d16 adds pairs 8 and 9, which
  bring one new daughter each.
- The largest network (d16) has 17 modes and 9 triplets.

**Outcome:** every N, with or without closure, gives the same run:

| ⟨E_parents⟩ | swing | dissipation/driving | active daughters | stop |
|---|---|---|---|---|
| 1.498e-10 | 0.008 | 1.0005 | 2, (3,−31) and (5,−8) | end, t = 10¹³ s |

The parent settles at the E_eq of pair 1 (1.50e-10). That is below the E_th of every later pair (lowest:
1.30e-09), so pairs 2–9 never grow. Up to d16, the ladder only ever tests its first pair.

## `ladder_2p_shared` (and `ladder_2p_shared_direct`)

Parents (2,−3), f = 9.373 c/d, γ = −8.642e-13 s⁻¹, and (2,−2) as above.

| step | shared daughter d | pair under (2,−3) (E_th) | pair under (2,−2) (E_th) |
|---|---|---|---|
| 1 | (6,−12) | (6,−55) + (6,−12) (3.61e-07) | (6,−29) + (6,−12) (1.76e-07) |
| 2 | (3,−6) | (1,−16) + (3,−6) (1.33e-06) | (3,−6) + (5,−29) (6.42e-07) |
| 3 | (3,−10) | (3,−12) + (3,−10) (8.86e-08) | (3,−10) + (5,−14) (7.78e-06) |
| 4 | (6,−13) | (6,−13) + (8,−62) (4.41e-06) | (4,−18) + (6,−13) (7.79e-06) |
| 5 | (8,−35) | (6,−17) + (8,−35) (1.66e-08) | (6,−13) + (8,−35) (9.06e-06) |
| 6 | (3,−23) | (3,−23) + (5,−11) (1.18e-05) | (3,−23) + (5,−9) (6.23e-06) |
| 7 | (5,−11) | already in from step 6 | (5,−21) + (5,−11) (7.17e-08) |
| 8 | (6,−14) | (4,−26) + (6,−14) (4.93e-08) | (6,−14) + (8,−31) (1.19e-05) |

Daughters (f [c/d], γ [s⁻¹]):

| | | | |
|---|---|---|---|
| (6,−55) 1.897, 4.83e-09 | (6,−12) 7.477, 2.30e-12 | (6,−29) 3.470, 7.94e-11 | (1,−16) 1.359, 4.24e-12 |
| (3,−6) 8.009, 5.65e-13 | (5,−29) 2.952, 8.27e-11 | (3,−12) 4.230, 2.19e-12 | (3,−10) 5.144, 1.35e-12 |
| (5,−14) 5.791, 3.63e-12 | (6,−13) 7.156, 2.69e-12 | (8,−62) 2.208, 1.08e-08 | (4,−18) 3.807, 5.59e-12 |
| (6,−17) 5.597, 5.56e-12 | (8,−35) 3.776, 1.30e-10 | (3,−23) 2.367, 1.80e-11 | (5,−11) 6.992, 2.05e-12 |
| (5,−9) 8.548, 1.11e-12 | (5,−21) 3.952, 1.41e-11 | (4,−26) 2.671, 3.76e-11 | (6,−14) 6.702, 3.94e-12 |
| (8,−31) 4.235, 7.04e-11 | | | |

**Outcome:**

| network | modes | ⟨E_parents⟩ | swing | dissipation/driving | active daughters | stop |
|---|---|---|---|---|---|---|
| shared1_d3 | 5 | 5.38e-06 | 0.000 | 0.9997 | 2 | end |
| shared2_d6 | 8 | 1.25e-05 | 4.3 | 0.82 | 4 | end |
| shared3_d9 | 11 | 1.95e-04 | 9.7 | 0.28 | 5 | runaway, t = 1.06e12 s |
| shared4_d12 | 14 | 1.67e-04 | 20.8 | 0.11 | 5 | runaway, t = 1.01e12 s |
| shared5_d14 | 16 | 1.21e-04 | 26.5 | 0.15 | 6 | runaway, t = 1.50e12 s |
| shared6_d17 | 19 | 1.18e-04 | 29.7 | 0.15 | 6 | runaway, t = 1.50e12 s |
| shared7_d18 | 20 | 1.45e-06 | 4.7 | 0.96 | 9 | end, not settled |
| shared8_d21 | 23 | 1.41e-06 | 2.4 | 0.986 | 17 | end, not settled |

Only shared daughters enter this set, so (2,−2)'s own best pair, (3,−31) + (5,−8) with E_th = 4.9e-12, is
never in it. The parents therefore climb to 10⁻⁶ – 10⁻⁴ E★, 10⁴ – 10⁶ times their energy in `ladder_1p`,
before anything holds them.

## Updates (2026-11-11)

1. **Folder move.** `ladder_2p_direct` moved here from `out/networks_mixed/`, and future mixed runs write
   here. The other models' mixed sets moved to their own `out/networks/<model>/` folders.
2. **The `_direct` folders hold no direct coupling.**
   - Neither (2,−2) with itself nor (2,−2) with (2,−3) has a damped c within the cut.
   - So `ladder_1p_direct` and `ladder_2p_shared_direct` contain exactly the networks and runs of
     `ladder_1p` and `ladder_2p_shared`.
   - They are left over from builds before the builder stopped writing a `_direct` file without a direct
     triplet. They can be deleted.
3. **Closure changes nothing here.** `ladder_1p` with and without closure are the same networks: the same
   triplets, only listed in a different order.
4. **Direct ranking.**
   - These files were ranked by R_c; the builder now ranks by E_eq, and μ gives the same order as R_c.
   - All three pick the same c for both direct triplets above, so rebuilding gives the same networks; only
     the header columns change.
   - Across all 92 parent pairs of this model, E_eq and μ pick the same c.
5. **Plots redrawn.** The bottom panel is now labelled "Driving / damping power, 2|γ|E [E★ s⁻¹]".

## Open before the next step

- The mixed set never reaches a two-parent state: (1,+7)'s hold energy is about 50 × K_s. A mixed state
  needs a fast partner held below K_s of its slow parent. No candidate of this model has been checked
  against that yet.
- `ladder_1p` beyond d2 adds nothing at this parent, and `ladder_2p_shared` is missing (2,−2)'s best pair.
  Whether to change the ladder rule is open.
