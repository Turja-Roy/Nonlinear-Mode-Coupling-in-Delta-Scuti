#!/usr/bin/env bash
# Mixed networks (Plans/network-growth.md, "Sum-branch parents"): two driven
# parents x, y with a damped c near f_x + f_y, plus daughter pairs.
# One network_sweep job per model, sets 2p_direct (each parent its own pairs)
# and 2p_shared_direct (shared daughters, MW25 Fig. 6 type).
#
#   LIST_DIRECT=1 SETS=" " sbatch shell-scripts/network_sweep.sh   # 1. rank pairs, out/direct_candidates_*.txt
#   BUILD_ONLY=1 shell-scripts/mixed_submit.sh                      # 2. build, read the logs
#   shell-scripts/mixed_submit.sh                                   # 3. build + run
#   shell-scripts/mixed_submit.sh --time=4:00:00                    # extra args go to sbatch
#
# A model goes ahead only if its build log shows "1 direct" and each parent
# "N with E_eq" > 0. PAIRS: from the observables tables (gamma_rad, l <= 3,
# |n| <= 20), one current parent kept; replace with rows of step 1.
# Output: data/networks_mixed/, out/networks_mixed/.

set -euo pipefail
source shell-scripts/config.sh
declare -A PAIRS=(
  [dsct_M2.0]="2,-2 2,5"
  [dsct_M2.0_T7696]="2,-2 1,7"
  [dsct_M1.85_T7555]="1,-1 2,4"
  [dsct_M2.2_T7888]="1,3 1,4"
  [dsct_M1.7_T7750]="0,7 2,-3"
)
for i in "${!MODELS[@]}"; do
  tag=${MODELS[$i]%%:*}
  [[ -n ${PAIRS[$tag]:-} ]] || continue
  echo "$tag: parents ${PAIRS[$tag]}"
  PARENTS="${PAIRS[$tag]}" SETS="${SETS:-2p_direct 2p_shared_direct}" \
    DATA=data/networks_mixed OUT=out/networks_mixed \
    sbatch --array="$i" --job-name="mixed_${tag#dsct_}" "$@" shell-scripts/network_sweep.sh
done
