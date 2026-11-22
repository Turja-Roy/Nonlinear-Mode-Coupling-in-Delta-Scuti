#!/usr/bin/env bash
# Network ladders on the OSPool (ap40), HTCondor. Run it on the access point itself, from the repository:
#
#   SETS=2p_sum shell-scripts/osg/network_sweep.sh             # MODEL=models/<tag> for another model
#   BUILD_ONLY=1 SETS=1p shell-scripts/osg/network_sweep.sh    # build and place, no runs
#
# Sets, settings and outputs: shell-scripts/network_lib.sh. Compiles and builds on the access point
# (libraries: shell-scripts/osg/setup.sh, once), then submits one 1-core job per network still to run
# (shell-scripts/osg/network_job.sh; network_run with the shared libraries it links, glibc aside, in
# build/lib), so every network runs at once on as many machines as the pool gives. A DAG: its POST step,
# on the access point, plots the runs, then the convergence and networks.md (shell-scripts/network_post.sh,
# <out>/<model>/post.log). WALL [min] a run, default 1800: 30 h, below the 40 h of a "Long" job.
# Do not start it again while its DAG runs (condor_q): the run folders of the queued networks are emptied.

set -euo pipefail
source shell-scripts/config.sh
S=/ospool/ap40/data/$USER/soft
export SUNDIALS_DIR=$S/sundials HDF5_DIR=$S/mesasdk PATH=$S/gsl/bin:$PATH \
    PKG_CONFIG_PATH=$S/gsl/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}
setup_python
mkdir -p build
flock build/.compile.lock ./compile.sh -j 4 network_build network_run
! ldd build/network_run | grep 'not found' || exit 1
rm -rf build/lib; mkdir build/lib
ldd build/network_run | awk '$3 ~ /^\// {print $3}' | grep -vE '/(libc|libm|libpthread|libdl|librt|ld-linux[^/]*)\.so' \
    | xargs -r cp -L -t build/lib
JOBS=4
WALL=${WALL:-1800}
source shell-scripts/network_lib.sh
echo "=== $MODEL"
[[ -z ${SETS// } ]] && exit 0
build_sets
[[ ${BUILD_ONLY:-0} == 1 ]] && exit 0

J=$PWD/$NET_DATA/$TAG/condor
rm -rf "$J"; mkdir -p "$J"
todo | while read -r d; do realpath "$d"; done > "$J/jobs.txt"
if [[ ! -s $J/jobs.txt ]]; then shell-scripts/network_post.sh "$PWD" "$OUT/$TAG"; exit 0; fi
cat > "$J/runs.sub" << EOF
executable              = $PWD/shell-scripts/osg/network_job.sh
arguments               = "$WALL $RUN_ARGS"
initialdir              = \$(dir)
transfer_input_files    = network.data, $PWD/build/network_run, $PWD/build/lib
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
output                  = condor.out
error                   = condor.err
log                     = $J/runs.log
request_cpus            = 1
request_memory          = 1GB
request_disk            = 1GB
# the access point's OS, for its glibc
requirements            = OSGVO_OS_STRING == "RHEL 9"
+JobDurationCategory    = "Long"
queue dir from $J/jobs.txt
EOF
printf 'JOB runs runs.sub\nSCRIPT POST runs %s %s %s %s\n' "$PWD/shell-scripts/network_post.sh" "$PWD" \
    "$PWD/$OUT/$TAG" "$PWD/$OUT/$TAG/post.log" > "$J/sweep.dag"
echo "$(wc -l < "$J/jobs.txt") networks to run"
cd "$J" && condor_submit_dag sweep.dag
