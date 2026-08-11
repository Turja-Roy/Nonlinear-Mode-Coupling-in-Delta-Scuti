#!/usr/bin/env bash
# Shared configuration for the five-model run. Sourced by the job scripts.
#
# Conventions follow damsa-simulation/shell_scripts/submit_library.sh, which
# is known to submit successfully on this cluster: partition `normal`, no
# --account / --qos, miniforge + a repo-local .venv for Python, logs/ at the
# repo root, WORKDIR from SLURM_SUBMIT_DIR.

# ---------------------------------------------------------------- software
# Already installed on the cluster -- nothing is built by these jobs.
export MESASDK_ROOT="${MESASDK_ROOT:-$HOME/soft/mesasdk}"
export MESA_DIR="${MESA_DIR:-$HOME/soft/mesa-26.04.1}"
export GYRE_DIR="${GYRE_DIR:-$HOME/soft/gyre-9.1.1}"

# ---------------------------------------------------------------- pipeline
# Regenerate every model's GYRE inlists with scripts/make_inlists.py so all
# five are produced identically. This rewrites dsct_M2.0's hand-written
# inlists; the job copies them to *.in.handwritten first.
REGENERATE_INLISTS="${REGENERATE_INLISTS:-1}"

# Run the l = 16-25 second pass of the wide net. Set to 0 to stop after pass 1
# and inspect four-mode candidates first.
RUN_WIDE_PASS2="${RUN_WIDE_PASS2:-1}"

# The five models: TAG:MASS:LOGG:MW23_TEFF
# log g is what the model is selected on; Teff is reported, never asserted
# (Plans/nonlinear-coupling-pipeline.md sec 1.1).
MODELS=(
  "dsct_M2.0:2.0:3.900:7696"
  "dsct_M2.0_T7696:2.0:3.930:7696"
  "dsct_M2.2_T7888:2.2:3.830:7888"
  "dsct_M1.85_T7555:1.85:4.010:7555"
  "dsct_M1.7_T7750:1.7:4.230:7750"
)

model_field() {  # model_field <tag> <1=mass|2=logg|3=teff>
  local tag="$1" idx="$2" entry
  for entry in "${MODELS[@]}"; do
    if [ "${entry%%:*}" = "$tag" ]; then
      echo "$entry" | cut -d: -f$((idx + 1))
      return 0
    fi
  done
  echo "unknown model tag: $tag" >&2
  return 1
}

# ------------------------------------------------------------------ python
# Same pattern as damsa's submit_library.sh. The venv is rebuilt if it is
# missing or was copied from another machine (its interpreter will not run),
# which is what happens when the local .venv/ is rsynced up by accident.
setup_python() {
  unset PYTHONPATH
  unset PYTHONHOME
  # shellcheck disable=SC1090
  source ~/miniforge/bin/activate
  # All models share one .venv, so two jobs starting together could otherwise
  # rm -rf it out from under each other. Serialize the whole setup.
  exec 9>.venv.lock
  flock 9 2>/dev/null || echo "    (no flock; .venv setup unserialized)"
  if [ ! -x .venv/bin/python ] || ! .venv/bin/python -c '' 2>/dev/null; then
    echo "    (re)creating .venv"
    rm -rf .venv
    python -m venv .venv
  fi
  # shellcheck disable=SC1091
  source .venv/bin/activate
  # preflight is the real gate; a transient index failure must not kill a job
  # whose deps are already installed.
  pip install -q -r requirements.txt || echo "    pip install failed; continuing"
  exec 9>&-
  PYTHON="$(command -v python)"
  export PYTHON
}

setup_mesa_gyre() {
  if [ -f "$MESASDK_ROOT/bin/mesasdk_init.sh" ]; then
    # shellcheck disable=SC1091
    source "$MESASDK_ROOT/bin/mesasdk_init.sh"
  else
    echo "WARNING: no mesasdk at $MESASDK_ROOT" >&2
  fi
  export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
  ulimit -s unlimited 2>/dev/null || true   # MESA wants a big stack
}

preflight() {
  local fail=0
  [ -x "$GYRE_DIR/bin/gyre" ] || { echo "no gyre at $GYRE_DIR/bin/gyre" >&2; fail=1; }
  [ -d "$MESA_DIR/star" ]     || { echo "no MESA at $MESA_DIR" >&2; fail=1; }
  "$PYTHON" -c "import numpy, scipy, pandas, h5py, matplotlib" \
    || { echo "python deps missing" >&2; fail=1; }
  return $fail
}
