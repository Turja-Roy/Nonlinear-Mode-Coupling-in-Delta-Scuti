#!/bin/bash
#SBATCH -J dsct
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 48:00:00
#SBATCH -o logs/dsct_%j.out
#SBATCH -e logs/dsct_%j.err
#
# The whole pipeline for ONE model, end to end:
#
#   MESA -> select profile -> GYRE inlists -> GYRE narrow -> coupling tables
#        -> GYRE wide pass 1 -> pass 2 -> four-mode -> M6 -> figures
#
#   sbatch -J dsct_M2.2 shell-scripts/10_pipeline.sbatch dsct_M2.2_T7888
#
# submit_all.sh overrides -J and adds --dependency where needed. Every stage
# is skipped when its output already exists, so a job killed by walltime can
# be resubmitted verbatim and continues. Pass --force-all to redo everything.
#
# -c 16, no --mem: a wider request sits in PD Resources on this partition.
#   sacct -j <id> --format=JobName,Elapsed,MaxRSS,State

set -euo pipefail

WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$WORKDIR"
mkdir -p logs out

source shell-scripts/config.sh
setup_python
setup_mesa_gyre

TAG="${1:?usage: 10_pipeline.sbatch <model-tag> [--force-all]}"
FORCE_ALL="${2:-}"
MODEL="models/$TAG"
MASS=$(model_field "$TAG" 1)
LOGG=$(model_field "$TAG" 2)
TEFF=$(model_field "$TAG" 3)

echo "=== $TAG  (M=$MASS, log g=$LOGG, MW23 Teff=$TEFF) ==="
echo "    host=$(hostname)  cores=${SLURM_CPUS_PER_TASK:-$(nproc)}  started $(date)"
echo "    workdir=$WORKDIR"
echo "    MESA_DIR=$MESA_DIR"
echo "    GYRE_DIR=$GYRE_DIR"
test -d "$MODEL" || { echo "no such model dir: $MODEL" >&2; exit 2; }
preflight || exit 2

step() { echo; echo "--- [$TAG] $*"; }
have() { [ -s "$1" ] && [ -z "$FORCE_ALL" ]; }

# check_modes exits 2 on any failed check, and several of its tolerances were
# calibrated on dsct_M2.0 (2% nad miss rate, scan-edge margin). A different
# mass can trip one while being perfectly sound, so report and keep going --
# grep the log for FAILED before trusting a model's results.
check() {
  if ( cd "$MODEL/gyre" && "$PYTHON" check_modes.py "$@" ); then return 0; fi
  echo "!!! [$TAG] check_modes FAILED (${*:-narrow}) -- continuing anyway"
  FAILED_CHECKS="$FAILED_CHECKS ${*:-narrow}"
}
FAILED_CHECKS=""

# ---------------------------------------------------------------- 1. MESA
# A SHARES_TRACK marker means this model rides on another model's evolutionary
# track (dsct_M2.0_T7696 on dsct_M2.0's). submit_all.sh makes this job depend
# on that model's job, so the track is already there.
#
# The link is recreated here rather than shipped: the cluster mount does not
# preserve symlinks across rsync -- one written through it reads back empty
# and does not resolve -- so upload.sh skips it and the job makes its own.
# Everything downstream (load_model, make_inlists, select_model) expects a
# real models/<tag>/mesa path, hence a link rather than a variable.
if [ -f "$MODEL/SHARES_TRACK" ]; then
  OWNER=$(tr -d '[:space:]' < "$MODEL/SHARES_TRACK")
  test -d "models/$OWNER/mesa" || {
    echo "SHARES_TRACK names $OWNER but models/$OWNER/mesa does not exist" >&2
    exit 2
  }
  if [ ! -d "$MODEL/mesa" ]; then
    ln -sfn "../$OWNER/mesa" "$MODEL/mesa"
    echo "    linked $MODEL/mesa -> ../$OWNER/mesa"
  fi
  step "shares the $OWNER track; skipping MESA"
elif have "$MODEL/mesa/LOGS/history.data"; then
  step "MESA already run; skipping (--force-all to redo)"
else
  step "MESA"
  ( cd "$MODEL/mesa" && ./mk && ./rn )
fi

# ------------------------------------------------- 2. pick the target profile
step "selecting profile at log g = $LOGG"
SEL_OUT=$( cd "$MODEL/mesa" && "$PYTHON" select_model.py \
             --logg "$LOGG" --mass "$MASS" --teff "$TEFF" )
echo "$SEL_OUT"
PROFILE=$(echo "$SEL_OUT" | sed -n 's/^GYRE input: .*\/\(profile[0-9]*\.data\.GYRE\).*/\1/p')
test -n "$PROFILE" || { echo "could not parse selected profile" >&2; exit 3; }
# select_model.py still prints the path when the .GYRE twin is absent, marked
# "<-- MISSING"; do not hand GYRE a file that is not there.
test -s "$MODEL/mesa/LOGS/$PROFILE" || {
  echo "selected profile has no .GYRE twin: $MODEL/mesa/LOGS/$PROFILE" >&2
  echo "  (MESA needs write_pulse_data_with_profile = .true.)" >&2
  exit 3
}
echo "    -> $PROFILE"

# ------------------------------------------------------------- 3. GYRE inlists
if [ "$REGENERATE_INLISTS" = "1" ]; then
  step "generating GYRE inlists"
  for f in "$MODEL"/gyre/gyre_*.in; do
    if [ -e "$f" ] && [ ! -e "$f.handwritten" ]; then
      cp "$f" "$f.handwritten"
      echo "    backed up $(basename "$f") -> $(basename "$f").handwritten"
    fi
  done
  "$PYTHON" scripts/make_inlists.py "$MODEL" --profile "$PROFILE" --force
else
  step "using existing inlists"
fi

# Also done by make_inlists.py, but REGENERATE_INLISTS=0 skips that.
mkdir -p "$MODEL/gyre/detail" "$MODEL/gyre/detail_wide"

GYRE="$GYRE_DIR/bin/gyre"
run_gyre() {  # run_gyre <inlist> <summary-that-proves-it-ran>
  if have "$MODEL/gyre/$2"; then
    echo "    $1: $2 exists, skipping"
  else
    ( cd "$MODEL/gyre" && "$GYRE" "$1" )
  fi
}

# ------------------------------------------------------------ 4. narrow net
step "GYRE narrow net (l <= 3, |n_pg| <= 20)"
run_gyre gyre_ad.in  summary_ad.h5
run_gyre gyre_nad.in summary_nad.h5
check

# ------------------------------------------------------ 5. coupling tables
step "coupling tables: triplets, kappa, observables"
if have "out/observables_$TAG.csv"; then
  echo "    out/observables_$TAG.csv exists, skipping"
else
  "$PYTHON" scripts/three_mode_search.py --model "$MODEL" --force
fi

# ------------------------------------------------------- 6. wide net, pass 1
step "GYRE wide net pass 1 (l <= 15)"
run_gyre gyre_ad_wide.in  summary_ad_wide.h5
run_gyre gyre_nad_wide.in summary_nad_wide.h5
check --inlist gyre_ad_wide.in \
    --ad summary_ad_wide.h5 --nad summary_nad_wide.h5 --detail detail_wide
df -h "$MODEL/gyre" | tail -1

# ------------------------------------------------------- 7. wide net, pass 2
if [ "$RUN_WIDE_PASS2" = "1" ]; then
  step "GYRE wide net pass 2 (l = 16-25)"
  run_gyre gyre_ad_wide_hi.in  summary_ad_wide_hi.h5
  run_gyre gyre_nad_wide_hi.in summary_nad_wide_hi.h5
  check --inlist gyre_ad_wide_hi.in \
      --ad summary_ad_wide_hi.h5 --nad summary_nad_wide_hi.h5 --detail detail_wide
else
  step "wide net pass 2 disabled (RUN_WIDE_PASS2=0)"
fi

# ------------------------------------------------------------- 8. four-mode
step "four-mode search"
if have "out/four_mode_$TAG.csv"; then
  echo "    out/four_mode_$TAG.csv exists, skipping"
else
  "$PYTHON" scripts/four_mode_search.py --model "$MODEL" \
      --out "out/four_mode_$TAG.csv" -j "${SLURM_CPUS_PER_TASK:-$(nproc)}"
fi

# Both read CSVs that are already on disk by now, so a failure here costs a
# rerun of seconds, not of the wide net. Do not throw the job away over it.
step "M6 network integration"
"$PYTHON" scripts/m6_network.py --csv "out/four_mode_$TAG.csv" \
    --out "out/m6_network_$TAG.csv" \
  || { echo "!!! [$TAG] m6_network FAILED"; FAILED_CHECKS="$FAILED_CHECKS m6"; }

# --------------------------------------------------------------- 9. figures
step "figures"
"$PYTHON" scripts/make_figures.py --model "$MODEL" \
  || { echo "!!! [$TAG] make_figures FAILED"; FAILED_CHECKS="$FAILED_CHECKS figures"; }

echo
if [ -n "$FAILED_CHECKS" ]; then
  echo "=== $TAG done $(date) WITH FAILURES:$FAILED_CHECKS ==="
else
  echo "=== $TAG done $(date) ==="
fi
