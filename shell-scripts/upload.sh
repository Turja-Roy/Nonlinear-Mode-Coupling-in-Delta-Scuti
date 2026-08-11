#!/bin/bash
# Push the repo to the cluster mount, excluding everything the jobs regenerate.
#
#   bash shell-scripts/upload.sh                       # to the default mount
#   bash shell-scripts/upload.sh /some/other/path      # elsewhere
#   DRY=1 bash shell-scripts/upload.sh                 # list what would transfer
#
# The exclusions matter for more than transfer size: 10_pipeline.sbatch skips
# any stage whose output already exists, so shipping out/*.csv or a populated
# models/*/mesa/LOGS would make the dsct_M2.0 job skip nearly everything.
# Leaving them behind is what makes the cluster run a genuine rerun.
#
# gyre/detail{,_wide} are excluded for size; the job mkdirs them, since GYRE
# aborts rather than creating its own detail_template directory.
#
# Also excluded: .venv (its interpreter paths are absolute and would not run
# there -- config.sh recreates it), MESA build products (stale .o/.mod/.smod
# from a different toolchain break ./mk), restart_photo (written live by a
# running MESA job, and the local one names photos that are not there), and
# the Reading/ symlinks to a local books directory absent on the cluster.

set -euo pipefail
cd "$(dirname "$0")/.."

DEST="${1:-/home/turja/cluster_mount/Stellar}"

# -a is -rlptgoD; drop -g and -o. The cluster mount is a network filesystem
# that refuses chgrp/chown, so preserving group and owner fails on every path
# and rsync exits 23 even though all files transferred fine.
RSYNC_OPTS=(-rlptD --info=stats1,progress2 --human-readable)
[ "${DRY:-0}" = "1" ] && RSYNC_OPTS+=(--dry-run --itemize-changes)

mkdir -p "$DEST"

# Shared-track models keep a symlink at models/<tag>/mesa pointing at the model
# that owns the track. The mount does not preserve symlinks -- one written
# through it reads back empty and fails to resolve -- so do not send them.
# 10_pipeline.sbatch recreates the link from the SHARES_TRACK marker instead.
LINK_EXCLUDES=()
for marker in models/*/SHARES_TRACK; do
  [ -e "$marker" ] || continue
  LINK_EXCLUDES+=(--exclude="$(dirname "$marker")/mesa")
done
[ ${#LINK_EXCLUDES[@]} -gt 0 ] && \
  echo "skipping unshippable symlinks: ${LINK_EXCLUDES[*]//--exclude=/}"

rsync "${RSYNC_OPTS[@]}" \
  ${LINK_EXCLUDES[@]+"${LINK_EXCLUDES[@]}"} \
  --exclude='.git/' \
  --exclude='.venv/' \
  --exclude='__pycache__/' \
  --exclude='.pytest_cache/' \
  --exclude='.Old_Codes/' \
  --exclude='Reading/' \
  --exclude='logs/' \
  --exclude='out/' \
  --exclude='models/*/mesa/LOGS/' \
  --exclude='models/*/mesa/photos/' \
  --exclude='models/*/mesa/.mesa_temp_cache/' \
  --exclude='models/*/mesa/make/*.o' \
  --exclude='models/*/mesa/make/*.mod' \
  --exclude='models/*/mesa/make/*.smod' \
  --exclude='models/*/mesa/restart_photo' \
  --exclude='models/*/mesa/star' \
  --exclude='models/*/mesa/*.mod' \
  --exclude='models/*/gyre/detail/' \
  --exclude='models/*/gyre/detail_wide/' \
  --exclude='models/*/gyre/*.h5' \
  ./ "$DEST/"

echo
echo "uploaded to $DEST"
echo
echo "then, on the cluster:"
echo "    cd ~/Stellar          # or wherever $DEST maps to"
echo "    bash shell-scripts/submit_all.sh --dry-run"
echo "    bash shell-scripts/submit_all.sh"
