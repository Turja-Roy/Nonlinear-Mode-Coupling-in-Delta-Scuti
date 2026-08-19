#!/bin/bash
# Submit one pipeline job per model.
#
#   bash shell-scripts/submit_all.sh                       # all five
#   bash shell-scripts/submit_all.sh --dry-run             # print sbatch lines only
#   bash shell-scripts/submit_all.sh dsct_M2.2_T7888       # a subset
#
# MESA and GYRE are already installed under ~/soft, so there is no build job.
#
# A model carrying a SHARES_TRACK file rides on another model's evolutionary
# track, so its job depends on that model's job rather than running MESA twice.
# dsct_M2.0_T7696 shares dsct_M2.0's track this way.

set -euo pipefail
cd "$(dirname "$0")/.."          # repo root; becomes SLURM_SUBMIT_DIR
source shell-scripts/config.sh
mkdir -p logs out

DRY=0
WANTED=()
for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY=1 ;;
    -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
    *)         WANTED+=("$arg") ;;
  esac
done

submit() {  # submit <deps-or-empty> <sbatch-args...> ; echoes the job id
  local dep="$1"; shift
  # Every sbatch option must precede the script path -- anything after it is
  # passed to the script as a positional argument instead.
  local args=()
  # afterany, not afterok: the dependent only needs the owner's MESA track,
  # which is step 1 of 9. 10_pipeline.sbatch re-checks that the track is
  # really there, so a late failure in the owner must not cancel it.
  [ -n "$dep" ] && args+=(--dependency="afterany:$dep")
  args+=("$@")
  if [ "$DRY" = "1" ]; then
    echo "sbatch --parsable ${args[*]}" >&2
    echo "DRYRUN"
  else
    sbatch --parsable "${args[@]}"
  fi
}

declare -A JOB_OF
# Order matters: a model that shares a track must be submitted after its owner.
ORDER=()
SHARED=()
for entry in "${MODELS[@]}"; do
  tag="${entry%%:*}"
  if [ ${#WANTED[@]} -gt 0 ] && [[ ! " ${WANTED[*]} " =~ " $tag " ]]; then continue; fi
  if [ -f "models/$tag/SHARES_TRACK" ]; then SHARED+=("$tag"); else ORDER+=("$tag"); fi
done
[ ${#ORDER[@]} -eq 0 ] && [ ${#SHARED[@]} -eq 0 ] && { echo "no models selected" >&2; exit 2; }

for tag in ${ORDER[@]+"${ORDER[@]}"} ${SHARED[@]+"${SHARED[@]}"}; do
  dep=""
  if [ -f "models/$tag/SHARES_TRACK" ]; then
    owner=$(tr -d '[:space:]' < "models/$tag/SHARES_TRACK")
    if [ -s "models/$owner/mesa/LOGS/history.data" ]; then
      echo "  ($tag: $owner's track is already on disk -> no dependency)"
    elif [ -n "${JOB_OF[$owner]:-}" ]; then
      dep="${JOB_OF[$owner]}"
      echo "  ($tag shares the $owner track -> waits for job $dep)"
    else
      echo "  WARNING: $tag shares $owner's track but $owner was not submitted;" >&2
      echo "           make sure models/$owner/mesa/LOGS is already populated." >&2
    fi
  fi
  id=$(submit "$dep" -J "$tag" -o "logs/${tag}_%j.out" -e "logs/${tag}_%j.err" \
              shell-scripts/10_pipeline.sbatch "$tag")
  JOB_OF[$tag]="$id"
  printf '%-18s -> job %s\n' "$tag" "$id"
done

echo
if [ "$DRY" = "1" ]; then
  echo "dry run: nothing submitted"
else
  echo "watch : squeue -u \$USER"
  echo "logs  : logs/<tag>_<jobid>.{out,err}"
fi
