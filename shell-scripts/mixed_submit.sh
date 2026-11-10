#!/usr/bin/env bash
# Mixed networks (Plans/network-growth.md, "Sum-branch parents"): two driven
# parents x, y with a damped c near f_x + f_y, plus daughter pairs.
# One network_sweep job per model, set 2p_direct (each parent its own pairs).
# 2p_shared_direct (MW25 Fig. 6 type) builds nothing here: the fast parent has
# no damped pairs, so no daughter is shared; SETS="2p_direct 2p_shared_direct" to try.
#
#   LIST_DIRECT=1 SETS=" " sbatch shell-scripts/network_sweep.sh   # 1. rank pairs, out/direct_candidates_*.txt
#   BUILD_ONLY=1 shell-scripts/mixed_submit.sh                      # 2. build, read the logs
#   shell-scripts/mixed_submit.sh                                   # 3. build + run
#   shell-scripts/mixed_submit.sh --time=4:00:00                    # extra args go to sbatch
#
# PAIRS, from out/direct_candidates_*.txt (2026-11-04): no x + y -> c in any model
# has both parents with an E_eq pair (two slow modes sum into the driven band).
# So: slow parent (E_eq pairs) + fast high-order p mode (no damped pairs), the
# fast one held by its harmonic x + x -> c' (zeta = |gamma_x|/gamma_c' < 1/2)
# and by x + y -> c once E_slow is large enough. M2.2: no harmonic hold; its
# fast (1,4) is held by the sum once E_(1,3) > 4.4e-10 (its E_th 3.6e-10).
# Output: data/networks_mixed/, out/networks_mixed/.

set -euo pipefail
source shell-scripts/config.sh
declare -A PAIRS=(
  [dsct_M2.0]="2,-2 2,5"
  [dsct_M2.0_T7696]="2,-2 1,7"
  [dsct_M1.85_T7555]="0,1 2,6"
  [dsct_M2.2_T7888]="1,3 1,4"
  [dsct_M1.7_T7750]="0,1 2,6"
)
for i in "${!MODELS[@]}"; do
  tag=${MODELS[$i]%%:*}
  [[ -n ${PAIRS[$tag]:-} ]] || continue
  echo "$tag: parents ${PAIRS[$tag]}"
  PARENTS="${PAIRS[$tag]}" SETS="${SETS:-2p_direct}" \
    DATA=data/networks_mixed OUT=out/networks_mixed \
    sbatch --array="$i" --job-name="mixed_${tag#dsct_}" "$@" shell-scripts/network_sweep.sh
done
