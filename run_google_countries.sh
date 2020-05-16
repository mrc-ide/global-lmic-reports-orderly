#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
DATE=${1:-$TODAY}

echo "*** Date: $DATE"

echo "*** Updating country list"
./update_run_sh.R $DATE

echo "*** Running country reports"

# Parallel
grep -E '^[A-Z]{3}\s*' countries | \
parallel --progress -j 16 ./orderly run lmic_reports_google iso3c={} date=$DATE

