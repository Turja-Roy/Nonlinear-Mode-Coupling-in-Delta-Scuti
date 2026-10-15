#!/usr/bin/env bash
#SBATCH --job-name=network_sweep
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/network_sweep_%A_%a.log
# Network growth ladder (Plans/mode-networks.md, sec. 3-4), one array task per
# model in config.sh's MODELS, all five in parallel.
#
#   sbatch shell-scripts/network_sweep.sh            # parents: l <= 2 GYRE modes of lowest E_th
#   PARENTS="6,3 6,4" sbatch --array=2 shell-scripts/network_sweep.sh   # (l,n) differ per model
#   PARENT_LMAX=4 sbatch shell-scripts/network_sweep.sh   # let l <= 4 modes compete as parents
#   WALL=10 CLOSURE=0.15 sbatch shell-scripts/network_sweep.sh   # tighter cap, full closure
#   sbatch --array=2 shell-scripts/network_sweep.sh             # MODELS[2] only
#   NS="1 2 4 8 16" WALL=150 sbatch --cpus-per-task=32 shell-scripts/network_sweep.sh   # up to N16, 2.5 h cap
#
# Compiles network_build/network_run first (array tasks take turns on a lock;
# up-to-date objects are skipped), rebuilds every network from scratch (one
# model load per parent set), then integrates them in parallel, plotting each
# run as it finishes. A run is redone when its network file changed.
# Output: data/networks/<model>/<set>/, out/networks/<model>/<set>/<run>/.

set -euo pipefail
source shell-scripts/config.sh
setup_python
mkdir -p build
flock build/.compile.lock ./compile.sh -j "${SLURM_CPUS_PER_TASK:-4}" network_build network_run
MODEL=${MODEL:-models/$(cut -d: -f1 <<< "${MODELS[${SLURM_ARRAY_TASK_ID:-0}]}")}
echo "=== $MODEL (array task ${SLURM_ARRAY_TASK_ID:-none})"
TAG=$(basename "$MODEL")
JOBS=${SLURM_CPUS_PER_TASK:-4}
OUT=${OUT:-out/networks}
NS=${NS:-"1 2 4 8 16 32"}
CUTS="0.15 0.05 0.015 0.005 0.0015"
# Cost ~ t_end * max|Delta| * n_triplets / 4e5 s, so the closure cut on
# far-detuned couplings is the lever (selecting legs are always kept). t_end
# 1e13 s: a bounded 8-daughter net has settled by then. WALL [min]: a run
# stops itself at 90% of it (--max-seconds) and keeps what it has (stop =
# wall in summary.csv); `timeout` at WALL is only the backstop, logged as
# TIMEOUT with no summary. Keep WALL below the job's --time, or Slurm kills the
# job first and the plots never run.
# rtol 1e-8 matches 1e-10 to ~3e-4 of peak energies at 2-5x less cost. atol
# 1e-40 puts the numerical floor at E ~ 1e-80 (atol^2, every solver): regrowth
# after a deeper dip starts from that floor, so its timing depends on atol.
WALL=${WALL:-660}
CLOSURE=${CLOSURE:-0.005}
RUN_ARGS=${RUN_ARGS:-"--q0 1e-6 --q0-daughter 1e-9 --seed 1 --e-max 1e-2 --t-end 1e13 --rtol 1e-8 --atol 1e-40"}
export OMP_NUM_THREADS=1

# Parents: PARENTS="l,n l,n" if given, else the builder picks the driven l <= PARENT_LMAX (default 2)
# GYRE modes of lowest parametric threshold (one for the one-parent runs, two for the two-parent runs).
read -r P1 P2 <<< "${PARENTS:-}"
ONE=${P1:+--parent $P1}; ONE=${ONE:---n-parents 1}
TWO=${P2:+--parent $P1 --parent $P2}; TWO=${TWO:---n-parents 2}

rm -rf "data/networks/$TAG"          # stale files from older builds would be run too
build() {   # out-dir, build args...
    local d=$1; shift
    ./build/network_build --model "$MODEL" --kappa-cache "out/four_mode_${TAG}.kappa_cache.tsv" \
        -j "$JOBS" --out "data/networks/$TAG/$d" ${CLOSURE:+--closure-cut $CLOSURE} \
        ${PARENT_LMAX:+--parent-l-max $PARENT_LMAX} "$@" || echo "BUILD FAILED $d"
}   # e.g. --min-gamma-ratio leaves no pair; the other sets still run

# one parent: rules x N at the search cut, closure on and off; then the cut sweep
build one_parent $ONE --rule eth eeq delta dgamma random -N $NS
build one_parent $ONE --rule eth -N $NS --no-closure
build one_parent_cut $ONE --rule delta eth -N 16 --cut $CUTS
# two parents, taking turns for pairs
build two_parents $TWO --rule eth eeq delta -N $NS
build two_parents $TWO --rule eth -N $NS --no-closure
# drop pairs with a weakly damped daughter, min gamma_d >= r |gamma_p|
for r in 0.3 1; do
    build one_parent_g$r $ONE --rule eth eeq -N $NS --min-gamma-ratio $r
    build two_parents_g$r $TWO --rule eth eeq -N $NS --min-gamma-ratio $r
done
# each parent's best a -> d + d forced in
build two_parents_sc $TWO --rule eth eeq -N $NS --self-coupled 1

run() {   # network file
    local f=$1 d grp p=1
    grp=$(basename "$(dirname "$f")")
    d=$OUT/$TAG/$grp/$(basename "$f" .data)
    [[ $grp == two_parents* ]] && p=2
    [[ -f $d/summary.csv ]] && cmp -s "$f" "$d/network.data" && return 0
    rm -rf "$d"; mkdir -p "$d"; cp "$f" "$d/network.data"
    timeout "${WALL}m" ./build/network_run "$f" --out "$d" $RUN_ARGS --max-seconds $((WALL * 54)) \
        > "$d/log.txt" 2>&1
    case $? in
        0) python scripts/network_plot.py --data "$d" --parents $p --title "$TAG $(basename "$d")" \
               >> "$d/log.txt" 2>&1 || echo "PLOT FAILED $d" ;;
        124) echo "TIMEOUT $f (> ${WALL} min)" ;;
        *) echo "FAILED $f (see $d/log.txt)" ;;
    esac
}
export -f run
export OUT RUN_ARGS WALL TAG
# costliest first, ~ max|Delta| * n_triplets (Delta = col 5 * |omega_id1|), so no long run finishes alone
for f in data/networks/"$TAG"/*/*.data; do
    awk '/^#/ {next} NF==7 {w[$1]=$6} NF==5 {d=$5*w[$1]; if (d<0) d=-d; if (d>m) m=d; n++}
         END {print m*n, FILENAME}' "$f"
done | sort -gr | cut -d' ' -f2- | xargs -P "$JOBS" -I{} bash -c 'run {}'

for d in "$OUT/$TAG"/*/; do python scripts/network_convergence.py --runs "$d"; done
