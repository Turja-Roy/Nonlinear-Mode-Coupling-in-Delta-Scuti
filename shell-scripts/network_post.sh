#!/usr/bin/env bash
# The last step of shell-scripts/network_sweep.sh, for one model's runs: a plot for every finished run
# without one (on the OSPool, every run: the jobs do not plot), the convergence plot of every ladder with a
# finished run, then networks.md. On the OSPool it is the POST script of the sweep's DAG, and its output
# goes to <out>/<model>/post.log.
#   shell-scripts/network_post.sh <repo> <out>/<model> [post.log]
cd "$1"
[[ -n ${3:-} ]] && exec >> "$3" 2>&1
PY=${PY:-$([[ -x .venv/bin/python ]] && echo .venv/bin/python || echo python)}
m=$(basename "$2")
for d in "$2"/ladder_*/*/d*/; do
    d=${d%/}
    [[ -f $d/summary.csv && ! -f $d/network_run.png ]] || continue
    $PY scripts/network_plot.py --data "$d" --title "$m ${d#"$2/"}" >> "$d/log.txt" 2>&1 || echo "PLOT FAILED $d"
done
for d in "$2"/ladder_*/*/; do
    compgen -G "$d/d*/summary.csv" > /dev/null || continue   # no finished run, or the old layout
    $PY scripts/network_convergence.py --runs "$d" > /dev/null || echo "CONVERGENCE PLOT FAILED $d"
done
$PY scripts/network_index.py --model "$2" || echo "INDEX FAILED"
exit 0   # the DAG succeeds whatever single runs did
