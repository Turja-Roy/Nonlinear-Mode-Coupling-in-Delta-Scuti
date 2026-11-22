# Network ladders (Plans/merged-ladders.md): settings, building and placing the ladders, and the networks
# still to run. Shared by the UTA HPC sweep (shell-scripts/network_sweep.sh, Slurm) and the OSPool sweep
# (shell-scripts/osg/network_sweep.sh, HTCondor); no scheduler here. Source it from the repository root,
# after shell-scripts/config.sh, with WALL [min] and JOBS (threads of network_build) set.
#
# Sets (SETS, default "1p 1p_eq 2p"): 1p (a ladder per driven l <= 2 mode), 1p_eq (only pairs with an
# E_eq, --only-eeq), 2p (the NPAIRS best-coupled pairs of parents; PARENTS="2,-2 1,7" names one pair),
# 2p_sum (each SUM_PARENTS x with its sum-triplet partner, --sum-partner; runs go to ladder_2p).
# PARENT_LMAX: largest parent l, default 2. NMAX: largest daughter count. LIST_DIRECT=1: rank x + y -> c
# into out/direct_candidates_<model>.txt. Default model: config.sh's MODELS[1], dsct_M2.0_T7696, or the
# Slurm array index.
# Output: $NET_DATA/<model>/<set>/<parents>/d<N>.data, $OUT/<model>/ladder_{1p,2p}/<parents>/d<N>/, and
# $OUT/<model>/networks.md, $NET_DATA/<model>/<set>/candidates.csv. A ladder goes to its parents' folder,
# or, if that folder holds different networks from another rule (or --only-eeq), to the first
# <parents>_ladder<k> (k >= 2) free or of the same rule. A finished run is redone only when its network
# (modes and triplets) changed.

MODEL=${MODEL:-models/$(cut -d: -f1 <<< "${MODELS[${SLURM_ARRAY_TASK_ID:-1}]}")}
TAG=$(basename "$MODEL")
KAP="out/four_mode_${TAG}.kappa_cache.tsv"
OUT=${OUT:-out/networks}
NET_DATA=${NET_DATA:-data/networks}   # not DATA: the OSPool login sets that to the data area
SETS=${SETS:-1p 1p_eq 2p}
NMAX=${NMAX:-24}
NPAIRS=${NPAIRS:-4}
# parents that never settle in dsct_M2.0_T7696's one-parent ladders (Plans/direct-triplets-fixed-points.md §9)
SUM_PARENTS=${SUM_PARENTS:-"0,1 1,-1 0,2 0,3 1,1 1,2 1,3 2,0 2,1 0,5 0,6 0,7 1,4 1,5 1,6 2,3 2,4"}
# Cost ~ t_end * max|Delta| * n_triplets / 4e5 s, so the closure cut on far-detuned couplings is the lever
# (selecting triplets are always kept). t_end 1e13 s: a bounded 8-daughter net has settled by then. A run
# stops itself at 90% of WALL (--max-seconds) and keeps what it has (stop = wall in summary.csv).
# rtol 1e-8 matched 1e-10 to ~3e-4 of peak energies under rk8pd; recheck under Adams. atol stays at
# network_run's 1e-300: no floor in reach. Starting energies: network_run's defaults, parents |q| = 1e-6,
# daughters --daughter-ratio below.
CLOSURE=${CLOSURE:-0.005}
RUN_ARGS=${RUN_ARGS:-"--seed 1 --e-max 1e-2 --t-end 1e13 --rtol 1e-8"}
export OMP_NUM_THREADS=1

build() {   # set, build args...
    local d=$1; shift
    rm -rf "${NET_DATA:?}/$TAG/$d"     # stale files from older builds would be run too
    # the lock: two jobs appending new kappa to one cache file would interleave lines
    flock "$KAP.lock" ./build/network_build --model "$MODEL" --kappa-cache "$KAP" -j "$JOBS" --out "$NET_DATA/$TAG/$d" \
        ${CLOSURE:+--closure-cut $CLOSURE} ${PARENT_LMAX:+--parent-l-max $PARENT_LMAX} --n-max "$NMAX" "$@" \
        || echo "BUILD FAILED $d"
}
want() { [[ " $SETS " == *" $1 "* ]]; }
body() { grep -v '^#' "$1"; }
rule() { grep -m1 '^# network_build' "$1" | sed 's/ daughters [0-9]*//'; }
same() {   # built ladder, run folder: every built network already run there (extra run folders ignored)
    local f
    for f in "$1"/d*.data; do cmp -s <(body "$f") <(body "$2/$(basename "$f" .data)/network.data") || return 1; done
}
place() {  # set, run group
    local p n o k r
    for p in "$NET_DATA/$TAG/$1"/*/; do
        p=${p%/}; [[ -d $p ]] || continue    # no ladders (build failed)
        n=$(basename "$p"); k=1; o=$OUT/$TAG/$2/$n
        r=$(ls "$p"/d*.data | head -1)
        while r2=$(ls "$o"/d*/network.data 2>/dev/null | head -1) && [[ -n $r2 ]] \
              && [[ $(rule "$r") != $(rule "$r2") ]] && ! same "$p" "$o"; do
            k=$((k + 1)); o=$OUT/$TAG/$2/${n}_ladder$k
        done
        if [[ $k -gt 1 ]]; then mv "$p" "$NET_DATA/$TAG/$1/${n}_ladder$k"; echo "PLACED $1/$n -> $2/${n}_ladder$k"; fi
    done
}
build_sets() {   # the sets of SETS, built and placed
    if [[ ${LIST_DIRECT:-0} == 1 ]]; then
        ./build/network_build --model "$MODEL" --kappa-cache "$KAP" -j "$JOBS" --list-direct \
            ${PARENT_LMAX:+--parent-l-max $PARENT_LMAX} | tee "out/direct_candidates_${TAG}.txt"
    fi
    local two=--pairs\ $NPAIRS
    [[ -n ${PARENTS:-} ]] && two=$(printf -- '--parent %s ' $PARENTS)
    want 1p && build ladder_1p --each-parent
    want 1p_eq && build ladder_1p_eq --each-parent --only-eeq
    want 2p && build ladder_2p $two
    want 2p_sum && build ladder_2p_sum $(printf -- '--sum-partner %s ' $SUM_PARENTS)
    want 1p && place ladder_1p ladder_1p
    want 1p_eq && place ladder_1p_eq ladder_1p
    want 2p && place ladder_2p ladder_2p
    want 2p_sum && place ladder_2p_sum ladder_2p
    return 0
}
rundir() {   # network file -> its run folder
    local grp p
    grp=$(basename "$(dirname "$(dirname "$1")")"); p=$(basename "$(dirname "$1")")
    case $grp in ladder_1p_eq) grp=ladder_1p ;; ladder_2p_sum) grp=ladder_2p ;; esac
    echo "$OUT/$TAG/$grp/$p/$(basename "$1" .data)"
}
todo() {   # every network of SETS not yet run with this network: its run folder emptied, network.data in it,
           # its path printed; costliest first, ~ max|Delta| * n_triplets (Delta = col 5 * |omega_id1|)
    local f d
    shopt -s nullglob
    for f in $(for s in $SETS; do printf '%s ' "$NET_DATA/$TAG/ladder_$s"/*/*.data; done); do
        awk '/^#/ {next} NF==7 {w[$1]=$6} NF==5 {d=$5*w[$1]; if (d<0) d=-d; if (d>m) m=d; n++}
             END {print m*n, FILENAME}' "$f"
    done | sort -gr | cut -d' ' -f2- | while read -r f; do
        d=$(rundir "$f")
        if [[ -f $d/summary.csv ]] && cmp -s <(body "$f") <(body "$d/network.data"); then echo "KEPT $d" >&2; continue; fi
        rm -rf "$d"; mkdir -p "$d"; cp "$f" "$d/network.data"; echo "$d"
    done
}
