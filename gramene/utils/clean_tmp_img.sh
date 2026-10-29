#!/bin/bash
# Monthly cleanup of generated images; run from root's crontab at 00:00 on the 1st of each month.
# Removes image files older than 30 days, then any dated subfolders left empty.
set -eu
IMGDIR=${IMGDIR:-/usr/local/ensembl-live/tmp/img/tmp/temporary}
[ -d "$IMGDIR" ] || exit 0

find "$IMGDIR" -type f -mtime +30 -delete
find "$IMGDIR" -mindepth 1 -type d -empty -delete
