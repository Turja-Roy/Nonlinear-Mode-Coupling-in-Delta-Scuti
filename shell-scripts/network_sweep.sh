#!/usr/bin/env bash
#SBATCH --job-name=network_sweep
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=logs/network_sweep_%j.log
# Network growth ladder on one model (Plans/mode-networks.md, sec. 3-4).
#
#   sbatch scripts/network_sweep.sh                  # parents: l <= 2 GYRE modes of lowest E_th
#   PARENTS="6,3 6,4" sbatch scripts/network_sweep.sh
#   WALL=10 CLOSURE=0.15 sbatch scripts/network_sweep.sh    # tighter cap, full closure
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
# Cost ~ t_end * max|Delta| * n_triplets / 4e5 s, so the closure cut on
# far-detuned couplings is the lever (selecting legs are always kept). t_end
# 1e13 s: a bounded 8-daughter net has settled by then. WALL [min] kills a run
# from outside; it leaves no summary and is logged as TIMEOUT.
WALL=${WALL:-20}
CLOSURE=${CLOSURE:-0.005}
RUN_ARGS=${RUN_ARGS:-"--q0 1e-6 --q0-daughter 1e-9 --seed 1 --e-max 1e-2 --t-end 1e13"}
export OMP_NUM_THREADS=1

# Parents: PARENTS="l,n l,n" if given, else the builder picks the driven l <= 2
# GYRE modes of lowest parametric threshold (one for the one-parent runs, two for the two-parent runs).
read -r P1 P2 <<< "${PARENTS:-}"
ONE=${P1:+--parent $P1}; ONE=${ONE:---n-parents 1}
TWO=${P2:+--parent $P1 --parent $P2}; TWO=${TWO:---n-parents 2}

build() {   # out-dir, build args...
    local d=$1; shift
    ./build/network_build --model "$MODEL" --kappa-cache "out/four_mode_${TAG}.kappa_cache.tsv" \
        -j "$JOBS" --out "data/networks/$d" ${CLOSURE:+--closure-cut $CLOSURE} "$@"
}

# one parent: rules x N at the search cut, closure on and off; then the cut sweep
build one_parent $ONE --rule eth delta dgamma random -N $NS
build one_parent $ONE --rule eth -N $NS --no-closure
build one_parent_cut $ONE --rule delta eth -N 16 --cut $CUTS
# two parents, shared daughters compete with private ones under E_th
build two_parents $TWO --rule eth delta -N $NS
build two_parents $TWO --rule eth -N $NS --no-closure

run() {   # network file
    local f=$1 d
    d=$OUT/$(basename "$(dirname "$f")")/$(basename "$f" .data)
    [[ -f $d/summary.csv ]] && return 0
    mkdir -p "$d"
    timeout "${WALL}m" ./build/network_run "$f" --out "$d" $RUN_ARGS > "$d/log.txt" 2>&1
    case $? in
        0) ;;
        124) echo "TIMEOUT $f (> ${WALL} min)" ;;
        *) echo "FAILED $f (see $d/log.txt)" ;;
    esac
}
export -f run
export OUT RUN_ARGS WALL
ls data/networks/*/*.data | xargs -P "$JOBS" -I{} bash -c 'run {}'

for d in "$OUT"/*/; do python3 scripts/network_convergence.py --runs "$d"; done
