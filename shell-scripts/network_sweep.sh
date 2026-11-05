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
# Network growth ladders (Plans/network-growth.md), one array task per model
# in config.sh's MODELS, all five in parallel.
#
#   sbatch shell-scripts/network_sweep.sh            # parents: fields 5, 6 of each MODELS entry
#   PARENTS="6,3 6,4" sbatch --array=2 shell-scripts/network_sweep.sh   # override, MODELS[2] only
#   PARENTS=" " PARENT_LMAX=4 sbatch shell-scripts/network_sweep.sh     # builder picks, l <= 4
#   NS1="2 4" NS2="1 2" WALL=60 sbatch shell-scripts/network_sweep.sh   # short ladders, 1 h cap
#   LIST_DIRECT=1 SETS=" " sbatch shell-scripts/network_sweep.sh   # rank parent pairs x + y -> c only
#   SETS="2p_direct" BUILD_ONLY=1 ...                          # those sets, build, no runs
# Sets (SETS, default all but 2p_direct): 1p 1p_noclosure 1p_direct 2p_shared 2p_shared_direct
# 2p_direct. A *_direct set writes a file only where a direct triplet exists.
# Mixed networks: shell-scripts/mixed_submit.sh.
#
# Compiles network_build/network_run first (array tasks take turns on a lock;
# up-to-date objects are skipped), rebuilds every network from scratch (one
# model load per set), then integrates them in parallel, plotting each run as
# it finishes. A run is redone when its network file changed.
# Output: $DATA/<model>/<set>/, $OUT/<model>/<set>/<network>/, LIST_DIRECT:
# out/direct_candidates_<model>.txt.

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
SETS=${SETS:-1p 1p_noclosure 1p_direct 2p_shared 2p_shared_direct}
NS1=${NS1:-$(seq 2 2 16)}           # 1 parent: exactly 2, 4, ..., 16 daughters
NS2=${NS2:-$(seq 1 8)}              # 2 parents: 1, ..., 8 shared daughters
NS2P=${NS2P:-$(seq 2 2 16)}         # 2 parents, own pairs: exactly 2, 4, ..., 16 daughters
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

# Parents: PARENTS="l,n l,n", else fields 5, 6 of the MODELS entry, else the builder picks the
# driven l <= PARENT_LMAX (default 2) GYRE modes of lowest E_th (one, then two).
read -r P1 P2 <<< "${PARENTS:-$(model_field "$TAG" 4 2>/dev/null || true) $(model_field "$TAG" 5 2>/dev/null || true)}"
ONE=${P1:+--parent $P1}; ONE=${ONE:---n-parents 1}
TWO=${P2:+--parent $P1 --parent $P2}; TWO=${TWO:---n-parents 2}

if [[ ${LIST_DIRECT:-0} == 1 ]]; then
    ./build/network_build --model "$MODEL" --kappa-cache "$KAP" -j "$JOBS" --list-direct \
        ${PARENT_LMAX:+--parent-l-max $PARENT_LMAX} | tee "out/direct_candidates_${TAG}.txt"
fi
[[ -z ${SETS// } ]] && exit 0
rm -rf "$DATA/$TAG"                 # stale files from older builds would be run too
build() {   # out-dir, build args...
    local d=$1; shift
    ./build/network_build --model "$MODEL" --kappa-cache "$KAP" -j "$JOBS" \
        --out "$DATA/$TAG/$d" ${CLOSURE:+--closure-cut $CLOSURE} \
        ${PARENT_LMAX:+--parent-l-max $PARENT_LMAX} "$@" || echo "BUILD FAILED $d"
}
want() { [[ " $SETS " == *" $1 "* ]]; }
# 1 parent, one pair per step; closure on and off; + its best direct triplet a + a -> c
want 1p && build ladder_1p $ONE -N $NS1
want 1p_noclosure && build ladder_1p $ONE -N $NS1 --no-closure
want 1p_direct && build ladder_1p_direct $ONE -N $NS1 --direct 1
# 2 parents, one shared daughter (with its two pairs) per step; + the best direct triplet per pair of parents
want 2p_shared && build ladder_2p_shared $TWO --shared -N $NS2
want 2p_shared_direct && build ladder_2p_shared_direct $TWO --shared -N $NS2 --direct 1
# 2 parents, each its own pairs (best first, in turn) + the best direct triplet: MW25's mixed network
want 2p_direct && build ladder_2p_direct $TWO -N $NS2P --direct 1
[[ ${BUILD_ONLY:-0} == 1 ]] && exit 0

run() {   # network file
    local f=$1 d grp
    grp=$(basename "$(dirname "$f")")
    d=$OUT/$TAG/$grp/$(basename "$f" .data)
    [[ -f $d/summary.csv ]] && cmp -s "$f" "$d/network.data" && return 0
    rm -rf "$d"; mkdir -p "$d"; cp "$f" "$d/network.data"
    timeout "${WALL}m" ./build/network_run "$f" --out "$d" $RUN_ARGS \
        --max-seconds $((WALL * 54)) > "$d/log.txt" 2>&1
    case $? in
        0) python scripts/network_plot.py --data "$d" --title "$TAG $grp $(basename "$d")" \
               >> "$d/log.txt" 2>&1 || echo "PLOT FAILED $d" ;;
        124) echo "TIMEOUT $f (> ${WALL} min)" ;;
        *) echo "FAILED $f (see $d/log.txt)" ;;
    esac
}
export -f run
export OUT RUN_ARGS WALL TAG
# costliest first, ~ max|Delta| * n_triplets (Delta = col 5 * |omega_id1|), so no long run
# finishes alone
shopt -s nullglob
for f in "$DATA/$TAG"/*/*.data; do
    awk '/^#/ {next} NF==7 {w[$1]=$6} NF==5 {d=$5*w[$1]; if (d<0) d=-d; if (d>m) m=d; n++}
         END {print m*n, FILENAME}' "$f"
done | sort -gr | cut -d' ' -f2- | xargs -P "$JOBS" -I{} bash -c 'run "$@"' _ {}

for d in "$OUT/$TAG"/ladder_*/; do
    [[ -d $DATA/$TAG/$(basename "$d") ]] || continue python scripts/network_convergence.py --runs "$d" || echo "CONVERGENCE PLOT FAILED $d"; done
