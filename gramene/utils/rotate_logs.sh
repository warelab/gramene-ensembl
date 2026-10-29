#!/bin/bash
# Weekly rotation of the Ensembl 116 Apache logs; run from root's crontab at 01:00 every Monday.
# After a run, logs/ holds the current week (live files) plus last week (*.log.YYYYMMDD.gz) only.
set -eu
LOGDIR=${LOGDIR:-/usr/local/ensembl-live/logs}
STAMP=$(date +%Y%m%d)
cd "$LOGDIR"

# drop the previous rotation, i.e. logs from two weeks ago
rm -f -- *.log.*

for f in access.log error.log preload-errors.log; do
  [ -s "$f" ] || continue
  # copy-truncate: Apache keeps writing (O_APPEND) to the same open file, no restart needed
  cp -p -- "$f" "$f.$STAMP"
  : > "$f"
  gzip -f -- "$f.$STAMP"
done
