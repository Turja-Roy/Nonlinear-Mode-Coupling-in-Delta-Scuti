#!/usr/bin/env bash
#SBATCH --job-name=network_sweep
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=network_sweep_%j.log
# Network growth ladder on one model (Plans/mode-networks.md, sec. 3-4).
#
#   sbatch scripts/network_sweep.sh                  # parents from the four-mode parent table
#   PARENTS="6,3 6,4" sbatch scripts/network_sweep.sh
#
# Builds every network first (one model load per parent set), then integrates
# them in parallel. Run from the repo root after ./compile.sh.

set -euo pipefail
MODEL=${MODEL:-models/dsct_M2.0}
TAG=$(basename "$MODEL")
JOBS=${SLURM_CPUS_PER_TASK:-4}
OUT=${OUT:-out/networks}
NS="1 2 4 8 16 32"
CUTS="0.15 0.05 0.015 0.005 0.0015"
RUN_ARGS=${RUN_ARGS:-"--q0 1e-6 --q0-daughter 1e-9 --seed 1 --e-max 1e-2"}
export OMP_NUM_THREADS=1

# Best row of the parent branch: the pair MW25's mechanism has the best chance on.
if [[ -z ${PARENTS:-} ]]; then
    PARENTS=$(python3 -c "
import pandas as pd
r = pd.read_csv('out/four_mode_${TAG}_parent.csv').iloc[0]
print(f'{r.l_a},{r.n_a} {r.l_b},{r.n_b}')")
fi
read -r P1 P2 <<< "$PARENTS"
echo "parents: $P1 $P2"

build() {   # out-dir, build args...
    local d=$1; shift
    ./build/network_build --model "$MODEL" --kappa-cache "out/four_mode_${TAG}.kappa_cache.tsv" \
        -j "$JOBS" --out "data/networks/$d" "$@"
}

# one parent: rules x N at the search cut, closure on and off; then the cut sweep
build one_parent --parent "$P1" --rule eth delta dgamma random -N $NS
build one_parent --parent "$P1" --rule eth -N $NS --no-closure
build one_parent_cut --parent "$P1" --rule delta eth -N 16 --cut $CUTS
# two parents, shared daughters compete with private ones under E_th
build two_parents --parent "$P1" --parent "$P2" --rule eth delta -N $NS
build two_parents --parent "$P1" --parent "$P2" --rule eth -N $NS --no-closure

run() {   # network file
    local f=$1 d
    d=$OUT/$(basename "$(dirname "$f")")/$(basename "$f" .data)
    [[ -f $d/summary.csv ]] && return 0
    mkdir -p "$d"
    ./build/network_run "$f" --out "$d" $RUN_ARGS > "$d/log.txt" 2>&1 \
        || echo "FAILED $f (see $d/log.txt)"
}
export -f run
export OUT RUN_ARGS
ls data/networks/*/*.data | xargs -P "$JOBS" -I{} bash -c 'run {}'

for d in "$OUT"/*/; do python3 scripts/network_convergence.py --runs "$d"; done
