#!/bin/bash
# One network on an OSPool execute point: a job of shell-scripts/osg/network_sweep.sh. network_run on
# network.data, with the shared libraries it was linked against (lib/); HTCondor returns the new files to
# the run folder, and the sweep's POST step plots them on the access point.
#   network_job.sh <wall [min]> <network_run arguments...>
wall=$1
shift
LD_LIBRARY_PATH=$PWD/lib ./network_run network.data --out . "$@" --max-seconds $((wall * 54)) > log.txt 2>&1
