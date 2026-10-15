#!/usr/bin/env bash
#SBATCH --job-name=solver_compare
#SBATCH --cpus-per-task=12
#SBATCH --time=06:00:00
#SBATCH --output=logs/solver_compare_%j.log
# Same networks, five integrator settings (amplitude.hpp): rk8pd at the
# default atol (1e-14 max|A(0)| = 1e-20), at atol 1e-30 and at rtol 1e-8,
# CVODE Adams, CVODE BDF (dense, analytic Jacobian). Every run stops itself
# after MAX_S seconds of integration and still writes its outputs and solver
# counters (stop = wall), so slow settings report where their time went.
#
#   sbatch shell-scripts/solver_compare.sh
#   MAX_S=7200 MODEL=models/dsct_M2.2_T7888 sbatch shell-scripts/solver_compare.sh
#
# Output: out/solver_compare/<model>/<network>/<setting>/ and, at the end of
# the log, one row per run: time, steps, failed steps, rhs calls, Jacobians.

set -euo pipefail
source shell-scripts/config.sh
setup_python
mkdir -p build logs
flock build/.compile.lock ./compile.sh -j "${SLURM_CPUS_PER_TASK:-4}" network_build network_run
MODEL=${MODEL:-models/$(cut -d: -f1 <<< "${MODELS[0]}")}
TAG=$(basename "$MODEL")
JOBS=${SLURM_CPUS_PER_TASK:-4}
DATA=data/solver_compare/$TAG
OUT=out/solver_compare/$TAG
MAX_S=${MAX_S:-3600}
RUN_ARGS=${RUN_ARGS:-"--q0 1e-6 --q0-daughter 1e-9 --seed 1 --e-max 1e-2 --t-end 1e13"}
export OMP_NUM_THREADS=1

rm -rf "$DATA"
build() {   # out-dir, build args...
    local d=$1; shift
    ./build/network_build --model "$MODEL" --kappa-cache "out/four_mode_${TAG}.kappa_cache.tsv" \
        -j "$JOBS" --out "$DATA/$d" --closure-cut 0.005 "$@"
}
build one_parent --n-parents 1 --rule dgamma eth -N 1 8
build two_parents --n-parents 2 --rule eth -N 8

run() {   # network file, setting
    local f=$1 s=$2 d args
    d=$OUT/$(basename "$(dirname "$f")")_$(basename "$f" .data)/$s
    case $s in
        rk8pd)        args="--solver rk8pd" ;;
        rk8pd_atol30) args="--solver rk8pd --atol 1e-30" ;;
        rk8pd_rtol8)  args="--solver rk8pd --rtol 1e-8" ;;
        adams)        args="--solver adams" ;;
        bdf)          args="--solver bdf" ;;
    esac
    rm -rf "$d"; mkdir -p "$d"
    # timeout only as a backstop: --max-seconds stops the run cleanly first
    timeout "$((MAX_S / 60 + 15))m" ./build/network_run "$f" --out "$d" $RUN_ARGS $args \
        --max-seconds "$MAX_S" > "$d/log.txt" 2>&1 \
        && python scripts/network_plot.py --data "$d" --title "$TAG $(basename "$f" .data) $s" >> "$d/log.txt" 2>&1 \
        || echo "FAILED $d" >> "$d/log.txt"
}
export -f run
export OUT RUN_ARGS MAX_S TAG
for f in "$DATA"/*/*.data; do
    for s in rk8pd rk8pd_atol30 rk8pd_rtol8 adams bdf; do echo "$f $s"; done
done | xargs -P "$JOBS" -n 2 bash -c 'run "$0" "$1"'

python - "$OUT" <<'EOF'
import pathlib, sys
import pandas as pd
cols = ["stop", "t_end", "seconds", "steps", "failed_steps", "rhs_calls", "jac_evals", "nonlin_iters"]
rows = []
for d in sorted(pathlib.Path(sys.argv[1]).glob("*/*")):
    try:
        s = pd.read_csv(d / "summary.csv").iloc[0]
        e = pd.read_csv(d / "energies.csv").drop(columns="t").to_numpy().min()
        rows.append([d.parent.name, d.name] + [s[c] for c in cols] + [e])
    except (FileNotFoundError, KeyError):
        rows.append([d.parent.name, d.name] + ["-"] * len(cols) + ["-"])
pd.set_option("display.width", 250)
print(pd.DataFrame(rows, columns=["network", "setting"] + cols + ["lowest_E"]).to_string(index=False))
EOF
