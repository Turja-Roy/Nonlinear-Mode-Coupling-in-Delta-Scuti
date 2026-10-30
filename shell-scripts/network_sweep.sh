#!/usr/bin/env bash
#SBATCH --job-name=network_sweep
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/network_sweep_%A_%a.log
# Network growth ladders (Plans/network-growth.md), one array task per model
# in config.sh's MODELS, all five in parallel.
#
#   sbatch shell-scripts/network_sweep.sh            # parents: fields 5, 6 of each MODELS entry
#   PARENTS="6,3 6,4" sbatch --array=2 shell-scripts/network_sweep.sh   # override, MODELS[2] only
#   PARENTS=" " PARENT_LMAX=4 sbatch shell-scripts/network_sweep.sh     # builder picks, l <= 4
#   F_P=1.1 sbatch shell-scripts/network_sweep.sh    # f_p of the sets outside the f_p scan
#   NS1="2 4" NS2="1 2" WALL=60 sbatch shell-scripts/network_sweep.sh   # short ladders, 1 h cap
#
# Compiles network_build/network_run first (array tasks take turns on a lock;
# up-to-date objects are skipped), rebuilds every network from scratch (one
# model load per set), then integrates them in parallel, plotting each run as
# it finishes. A run is redone when its network file changed.
# Output: data/networks/<model>/<set>/, out/networks/<model>/<set>/<network>_fp<f_p>/.

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
NS1=${NS1:-$(seq 2 2 16)}           # 1 parent: exactly 2, 4, ..., 16 daughters
NS2=${NS2:-$(seq 1 8)}              # 2 parents: 1, ..., 8 shared daughters
F_PS=${F_PS:-"0.5 0.9 1.1 2"}       # parent starting energy / E_th, scanned on ladder_1p
F_P=${F_P:-0.9}                     # and for every other set
# Cost ~ t_end * max|Delta| * n_triplets / 4e5 s, so the closure cut on
# far-detuned couplings is the lever (selecting legs are always kept). t_end
# 1e13 s: a bounded 8-daughter net has settled by then. WALL [min]: a run
# stops itself at 90% of it (--max-seconds) and keeps what it has (stop =
# wall in summary.csv); `timeout` at WALL is only the backstop, logged as
# TIMEOUT with no summary. Keep WALL below the job's --time, or Slurm kills the
# job first and the plots never run.
# rtol 1e-8 matched 1e-10 to ~3e-4 of peak energies under rk8pd; recheck under
# Adams. atol stays at network_run's 1e-300: no floor in reach. Starting
# energies: network_run's f_p E_th, f_d = 1e-3 (f_p from F_PS / F_P below).
WALL=${WALL:-660}
CLOSURE=${CLOSURE:-0.005}
RUN_ARGS=${RUN_ARGS:-"--seed 1 --e-max 1e-2 --t-end 1e13 --rtol 1e-8"}
export OMP_NUM_THREADS=1

# Parents: PARENTS="l,n l,n", else fields 5, 6 of the MODELS entry, else the builder picks the
# driven l <= PARENT_LMAX (default 2) GYRE modes of lowest E_th (one, then two).
read -r P1 P2 <<< "${PARENTS:-$(model_field "$TAG" 4 2>/dev/null || true) $(model_field "$TAG" 5 2>/dev/null || true)}"
ONE=${P1:+--parent $P1}; ONE=${ONE:---n-parents 1}
TWO=${P2:+--parent $P1 --parent $P2}; TWO=${TWO:---n-parents 2}

rm -rf "data/networks/$TAG"          # stale files from older builds would be run too
build() {   # out-dir, build args...
    local d=$1; shift
    ./build/network_build --model "$MODEL" --kappa-cache "$KAP" -j "$JOBS" \
        --out "data/networks/$TAG/$d" ${CLOSURE:+--closure-cut $CLOSURE} \
        ${PARENT_LMAX:+--parent-l-max $PARENT_LMAX} "$@" || echo "BUILD FAILED $d"
}
# 1 parent, one pair per step; closure on and off; + its best direct leg a + a -> c
build ladder_1p $ONE -N $NS1
build ladder_1p $ONE -N $NS1 --no-closure
build ladder_1p_direct $ONE -N $NS1 --direct 1
# 2 parents, one shared daughter (with its two pairs) per step; + the best direct leg per pair of parents
build ladder_2p_shared $TWO --shared -N $NS2
build ladder_2p_shared_direct $TWO --shared -N $NS2 --direct 1

run() {   # network file, f_p
    local f=$1 fp=$2 d grp
    grp=$(basename "$(dirname "$f")")
    d=$OUT/$TAG/$grp/$(basename "$f" .data)_fp$fp
    [[ -f $d/summary.csv ]] && cmp -s "$f" "$d/network.data" && return 0
    rm -rf "$d"; mkdir -p "$d"; cp "$f" "$d/network.data"
    timeout "${WALL}m" ./build/network_run "$f" --out "$d" --eth-parent "$fp" $RUN_ARGS \
        --max-seconds $((WALL * 54)) > "$d/log.txt" 2>&1
    case $? in
        0) python scripts/network_plot.py --data "$d" --title "$TAG $grp $(basename "$d")" \
               >> "$d/log.txt" 2>&1 || echo "PLOT FAILED $d" ;;
        124) echo "TIMEOUT $f f_p $fp (> ${WALL} min)" ;;
        *) echo "FAILED $f f_p $fp (see $d/log.txt)" ;;
    esac
}
export -f run
export OUT RUN_ARGS WALL TAG
# costliest first, ~ max|Delta| * n_triplets (Delta = col 5 * |omega_id1|), so no long run
# finishes alone; ladder_1p once per f_p in F_PS
for f in data/networks/"$TAG"/*/*.data; do
    awk '/^#/ {next} NF==7 {w[$1]=$6} NF==5 {d=$5*w[$1]; if (d<0) d=-d; if (d>m) m=d; n++}
         END {print m*n, FILENAME}' "$f"
done | sort -gr | cut -d' ' -f2- | while read -r f; do
    if [[ $f == */ladder_1p/* ]]; then for p in $F_PS; do echo "$f $p"; done; else echo "$f $F_P"; fi
done | xargs -P "$JOBS" -L 1 bash -c 'run "$@"' _

for d in "$OUT/$TAG"/ladder_*/; do python scripts/network_convergence.py --runs "$d" || echo "CONVERGENCE PLOT FAILED $d"; done
