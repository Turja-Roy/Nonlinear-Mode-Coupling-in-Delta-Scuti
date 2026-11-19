#!/usr/bin/env bash
#SBATCH -J network_sweep
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=30000      # node has 31868 MB; the default per-CPU memory x 16 asks for more
#SBATCH -t 12:00:00
#SBATCH -a 0-4
#SBATCH -o logs/network_sweep_%A_%a.log
# Network ladders (Plans/merged-ladders.md), one array task per model in config.sh's MODELS.
#
#   SETS=1p sbatch --array=1 shell-scripts/network_sweep.sh     # a ladder per driven l <= 2 mode
#   SETS=2p sbatch --array=1 shell-scripts/network_sweep.sh     # the NPAIRS best-coupled pairs of parents
#   PARENTS="2,-2 1,7" SETS=2p sbatch --array=1 ...             # one named pair instead
#   SETS=2p_sum sbatch --array=1 shell-scripts/network_sweep.sh # each SUM_PARENTS x with its sum-triplet partner
#   LIST_DIRECT=1 SETS=" " sbatch ...                           # rank x + y -> c only
# Sets (SETS, default "1p 1p_eq 2p"): 1p, 1p_eq (only pairs with an E_eq, --only-eeq), 2p, 2p_sum
# (--sum-partner; runs go to ladder_2p).
# BUILD_ONLY=1: build, no runs. PARENT_LMAX: largest parent l, default 2. NMAX: largest daughter count.
#
# Compiles network_build/network_run first (array tasks take turns on a lock; up-to-date objects are
# skipped), rebuilds the sets asked for, then integrates their networks in parallel, plotting each
# run as it finishes; sets not in SETS are not run. A finished run is redone only when its network
# (modes and triplets) changed.
# Output: $DATA/<model>/<set>/<parents>/d<N>.data, $OUT/<model>/ladder_{1p,2p}/<parents>/d<N>/, and
# $OUT/<model>/networks.md, $DATA/<model>/<set>/candidates.csv. A ladder goes to its parents' folder,
# or, if that folder holds different networks from another rule (or --only-eeq), to the first
# <parents>_ladder<k> (k >= 2) free or of the same rule. LIST_DIRECT: out/direct_candidates_<model>.txt.

set -euo pipefail
source shell-scripts/config.sh
setup_python
mkdir -p build
flock build/.compile.lock ./compile.sh -j "${SLURM_CPUS_PER_TASK:-4}" network_build network_run
MODEL=${MODEL:-models/$(cut -d: -f1 <<< "${MODELS[${SLURM_ARRAY_TASK_ID:-0}]}")}
echo "=== $MODEL (array task ${SLURM_ARRAY_TASK_ID:-none})"
TAG=$(basename "$MODEL")
KAP="out/four_mode_${TAG}.kappa_cache.tsv"
JOBS=${SLURM_CPUS_PER_TASK:-4}
OUT=${OUT:-out/networks}
DATA=${DATA:-data/networks}
SETS=${SETS:-1p 1p_eq 2p}
NMAX=${NMAX:-24}
NPAIRS=${NPAIRS:-4}
# parents that never settle in dsct_M2.0_T7696's one-parent ladders (Plans/direct-triplets-fixed-points.md §9)
SUM_PARENTS=${SUM_PARENTS:-"0,1 1,-1 0,2 0,3 1,1 1,2 1,3 2,0 2,1 0,5 0,6 0,7 1,4 1,5 1,6 2,3 2,4"}
# Cost ~ t_end * max|Delta| * n_triplets / 4e5 s, so the closure cut on
# far-detuned couplings is the lever (selecting triplets are always kept). t_end
# 1e13 s: a bounded 8-daughter net has settled by then. WALL [min]: a run
# stops itself at 90% of it (--max-seconds) and keeps what it has (stop =
# wall in summary.csv); `timeout` at WALL is only the backstop, logged as
# TIMEOUT with no summary. Keep WALL below the job's --time, or Slurm kills the
# job first and the plots never run.
# rtol 1e-8 matched 1e-10 to ~3e-4 of peak energies under rk8pd; recheck under
# Adams. atol stays at network_run's 1e-300: no floor in reach. Starting
# energies: network_run's defaults, parents |q| = 1e-6, daughters --daughter-ratio below.
WALL=${WALL:-660}
CLOSURE=${CLOSURE:-0.005}
RUN_ARGS=${RUN_ARGS:-"--seed 1 --e-max 1e-2 --t-end 1e13 --rtol 1e-8"}
export OMP_NUM_THREADS=1

if [[ ${LIST_DIRECT:-0} == 1 ]]; then
    ./build/network_build --model "$MODEL" --kappa-cache "$KAP" -j "$JOBS" --list-direct \
        ${PARENT_LMAX:+--parent-l-max $PARENT_LMAX} | tee "out/direct_candidates_${TAG}.txt"
fi
[[ -z ${SETS// } ]] && exit 0
build() {   # set, build args...
    local d=$1; shift
    rm -rf "${DATA:?}/$TAG/$d"         # stale files from older builds would be run too
    # the lock: two jobs appending new kappa to one cache file would interleave lines
    flock "$KAP.lock" ./build/network_build --model "$MODEL" --kappa-cache "$KAP" -j "$JOBS" --out "$DATA/$TAG/$d" \
        ${CLOSURE:+--closure-cut $CLOSURE} ${PARENT_LMAX:+--parent-l-max $PARENT_LMAX} --n-max "$NMAX" "$@" \
        || echo "BUILD FAILED $d"
}
want() { [[ " $SETS " == *" $1 "* ]]; }
TWO=--pairs\ $NPAIRS
[[ -n ${PARENTS:-} ]] && TWO=$(printf -- '--parent %s ' $PARENTS)
want 1p && build ladder_1p --each-parent
want 1p_eq && build ladder_1p_eq --each-parent --only-eeq
want 2p && build ladder_2p $TWO
want 2p_sum && build ladder_2p_sum $(printf -- '--sum-partner %s ' $SUM_PARENTS)

body() { grep -v '^#' "$1"; }
rule() { grep -m1 '^# network_build' "$1" | sed 's/ daughters [0-9]*//'; }
same() {   # built ladder, run folder: every built network already run there (extra run folders ignored)
    local f
    for f in "$1"/d*.data; do cmp -s <(body "$f") <(body "$2/$(basename "$f" .data)/network.data") || return 1; done
}
place() {  # set, run group
    local p n o k r
    for p in "$DATA/$TAG/$1"/*/; do
        p=${p%/}; n=$(basename "$p"); k=1; o=$OUT/$TAG/$2/$n
        r=$(ls "$p"/d*.data | head -1)
        while r2=$(ls "$o"/d*/network.data 2>/dev/null | head -1) && [[ -n $r2 ]] \
              && [[ $(rule "$r") != $(rule "$r2") ]] && ! same "$p" "$o"; do
            k=$((k + 1)); o=$OUT/$TAG/$2/${n}_ladder$k
        done
        if [[ $k -gt 1 ]]; then mv "$p" "$DATA/$TAG/$1/${n}_ladder$k"; echo "PLACED $1/$n -> $2/${n}_ladder$k"; fi
    done
}
want 1p && place ladder_1p ladder_1p
want 1p_eq && place ladder_1p_eq ladder_1p
want 2p && place ladder_2p ladder_2p
want 2p_sum && place ladder_2p_sum ladder_2p
[[ ${BUILD_ONLY:-0} == 1 ]] && exit 0

run() {   # network file
    local f=$1 d grp p
    grp=$(basename "$(dirname "$(dirname "$f")")"); p=$(basename "$(dirname "$f")")
    case $grp in ladder_1p_eq) grp=ladder_1p ;; ladder_2p_sum) grp=ladder_2p ;; esac
    d=$OUT/$TAG/$grp/$p
    d=$d/$(basename "$f" .data)
    [[ -f $d/summary.csv ]] && cmp -s <(body "$f") <(body "$d/network.data") && { echo "KEPT $d (same network)"; return 0; }
    rm -rf "$d"; mkdir -p "$d"; cp "$f" "$d/network.data"
    timeout "${WALL}m" ./build/network_run "$f" --out "$d" $RUN_ARGS \
        --max-seconds $((WALL * 54)) > "$d/log.txt" 2>&1
    case $? in
        0) python scripts/network_plot.py --data "$d" --title "$TAG ${d#"$OUT/$TAG/"}" \
               >> "$d/log.txt" 2>&1 || echo "PLOT FAILED $d" ;;
        124) echo "TIMEOUT $f (> ${WALL} min)" ;;
        *) echo "FAILED $f (see $d/log.txt)" ;;
    esac
}
export -f run body
export OUT RUN_ARGS WALL TAG
# costliest first, ~ max|Delta| * n_triplets (Delta = col 5 * |omega_id1|), so no long run
# finishes alone
shopt -s nullglob
for f in $(for s in $SETS; do printf '%s ' "$DATA/$TAG/ladder_$s"/*/*.data; done); do   # this job's sets only
    awk '/^#/ {next} NF==7 {w[$1]=$6} NF==5 {d=$5*w[$1]; if (d<0) d=-d; if (d>m) m=d; n++}
         END {print m*n, FILENAME}' "$f"
done | sort -gr | cut -d' ' -f2- | xargs -P "$JOBS" -I{} bash -c 'run "$@"' _ {}

for d in "$OUT/$TAG"/ladder_*/*/; do
    compgen -G "$d/d*/summary.csv" > /dev/null || continue   # no finished run, or the old layout
    python scripts/network_convergence.py --runs "$d" > /dev/null || echo "CONVERGENCE PLOT FAILED $d"
done
python scripts/network_index.py --model "$OUT/$TAG" || echo "INDEX FAILED"
