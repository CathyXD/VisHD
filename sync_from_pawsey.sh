#!/bin/bash
# Run this ON the Peter Mac cluster to pull the VisHD directory from Pawsey.
# Usage:
#   ./sync_from_pawsey.sh            # sync
#   ./sync_from_pawsey.sh --dry-run  # preview changes only, no files transferred
#   ./sync_from_pawsey.sh --delete   # also remove local files no longer present on Pawsey

set -euo pipefail

REMOTE_USER="pawsey1172"
REMOTE_HOST="pawsey_data"
REMOTE_PATH="/scratch/pawsey1172/sweng/VisHD/"
LOCAL_PATH="${HOME}/VisHD/"

RSYNC_OPTS=(-avzPru
  --exclude "slurmout/"
  --exclude "raw/"
  --exclude ".git/")

for arg in "$@"; do
  case "$arg" in
    --dry-run) RSYNC_OPTS+=(--dry-run) ;;
    --delete)  RSYNC_OPTS+=(--delete) ;;
    *) echo "Unknown option: $arg" >&2; exit 1 ;;
  esac
done

mkdir -p "$LOCAL_PATH"
rsync "${RSYNC_OPTS[@]}" "${REMOTE_HOST}:${REMOTE_PATH}" "$LOCAL_PATH"
