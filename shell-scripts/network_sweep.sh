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
#   LIST_DIRECT=1 SETS=" " sbatch ...                           # rank x + y -> c only
# Sets (SETS, default "1p 1p_eq 2p"): 1p, 1p_eq (only pairs with an E_eq, --only-eeq), 2p.
# BUILD_ONLY=1: build, no runs. PARENT_LMAX: largest parent l, default 2. NMAX: largest daughter count.
#
# Compiles network_build/network_run first (array tasks take turns on a lock; up-to-date objects are
# skipped), rebuilds the sets asked for, then integrates their networks in parallel, plotting each
# run as it finishes. A finished run is redone only when its network (modes and triplets) changed.
# Output: $DATA/<model>/<set>/<parents>/d<N>.data, $OUT/<model>/ladder_{1p,2p}/<parents>/d<N>/, and
# $OUT/<model>/networks.md. A 1p_eq ladder, and a 2p ladder whose parents already have a ladder
# built by an earlier rule, go to <parents>_ladder2. LIST_DIRECT: out/direct_candidates_<model>.txt.

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
# a 2p ladder whose parents already have a ladder of an earlier rule (no "in the order they entered" line)
for p in "$DATA/$TAG"/ladder_2p/*/; do
    o=$OUT/$TAG/ladder_2p/$(basename "$p")
    [[ -d $o ]] && ! grep -qs 'in the order they entered' "$o"/*/network.data && mv "$p" "${p%/}_ladder2"
done
[[ ${BUILD_ONLY:-0} == 1 ]] && exit 0

body() { grep -v '^#' "$1"; }
run() {   # network file
    local f=$1 d grp p
    grp=$(basename "$(dirname "$(dirname "$f")")"); p=$(basename "$(dirname "$f")")
    case $grp in ladder_1p_eq) d=$OUT/$TAG/ladder_1p/${p}_ladder2 ;; *) d=$OUT/$TAG/$grp/$p ;; esac
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
for f in "$DATA/$TAG"/*/*/*.data; do
    awk '/^#/ {next} NF==7 {w[$1]=$6} NF==5 {d=$5*w[$1]; if (d<0) d=-d; if (d>m) m=d; n++}
         END {print m*n, FILENAME}' "$f"
done | sort -gr | cut -d' ' -f2- | xargs -P "$JOBS" -I{} bash -c 'run "$@"' _ {}

for d in "$OUT/$TAG"/ladder_*/*/; do
    python scripts/network_convergence.py --runs "$d" > /dev/null || echo "CONVERGENCE PLOT FAILED $d"
done
python scripts/network_index.py --model "$OUT/$TAG" || echo "INDEX FAILED"
