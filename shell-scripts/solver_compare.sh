#!/usr/bin/env bash
#SBATCH --job-name=solver_compare
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH --output=logs/solver_compare_%j.log
# Same networks, four integrators (amplitude.hpp): rk8pd at the default atol
# (1e-14 max|A(0)| = 1e-20) and at 1e-30, CVODE Adams, CVODE BDF. If rk8pd's
# low-energy floor and burst timing move with atol, they are numerical.
#
#   sbatch shell-scripts/solver_compare.sh
#   MODEL=models/dsct_M2.2_T7888 sbatch shell-scripts/solver_compare.sh
#
# Output: out/solver_compare/<model>/<network>/<solver>/ (run as network_sweep.sh)
# and a wall-time / lowest-energy table at the end of the log.

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
WALL=${WALL:-660}
RUN_ARGS=${RUN_ARGS:-"--q0 1e-6 --q0-daughter 1e-9 --seed 1 --e-max 1e-2 --t-end 1e13"}
export OMP_NUM_THREADS=1

rm -rf "$DATA"
build() {   # out-dir, build args...
    local d=$1; shift
    ./build/network_build --model "$MODEL" --kappa-cache "out/four_mode_${TAG}.kappa_cache.tsv" \
        -j "$JOBS" --out "$DATA/$d" --closure-cut 0.005 "$@"
}
build one_parent --n-parents 1 --rule dgamma eth -N 1 32
build two_parents --n-parents 2 --rule eth -N 32

run() {   # network file, solver label
    local f=$1 s=$2 d args
    d=$OUT/$(basename "$(dirname "$f")")_$(basename "$f" .data)/$s
    case $s in
        rk8pd)       args="--solver rk8pd" ;;
        rk8pd_atol30) args="--solver rk8pd --atol 1e-30" ;;
        adams)       args="--solver adams" ;;
        bdf)         args="--solver bdf" ;;
    esac
    rm -rf "$d"; mkdir -p "$d"
    local t0=$SECONDS
    timeout "${WALL}m" ./build/network_run "$f" --out "$d" $RUN_ARGS $args > "$d/log.txt" 2>&1 \
        && python scripts/network_plot.py --data "$d" --title "$TAG $(basename "$f" .data) $s" >> "$d/log.txt" 2>&1 \
        || echo "FAILED or TIMEOUT $d" >> "$d/log.txt"
    echo "$((SECONDS - t0))" > "$d/wall_s.txt"
}
export -f run
export OUT RUN_ARGS WALL TAG
for f in "$DATA"/*/*.data; do
    for s in rk8pd rk8pd_atol30 adams bdf; do echo "$f $s"; done
done | xargs -P "$JOBS" -n 2 bash -c 'run "$0" "$1"'

python - "$OUT" <<'EOF'
import pathlib, sys
import pandas as pd
print(f"{'network':40s} {'solver':13s} {'wall [s]':>9s} {'rows':>5s}  lowest E of any mode")
for d in sorted(pathlib.Path(sys.argv[1]).glob("*/*")):
    wall = (d / "wall_s.txt").read_text().strip() if (d / "wall_s.txt").exists() else "?"
    try:
        df = pd.read_csv(d / "energies.csv")
        low = f"{df.drop(columns='t').to_numpy().min():.1e}"
        rows = len(df)
    except FileNotFoundError:
        low, rows = "no output", 0
    print(f"{d.parent.name:40s} {d.name:13s} {wall:>9s} {rows:5d}  {low}")
EOF
