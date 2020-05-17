#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
COUNTRY=$1
DATE=${2:-$TODAY}
SHORT_RUN=$3

echo "*** Country: $COUNTRY"
echo "*** Date: $DATE"

# Parallel
./orderly run lmic_reports_google iso3c=$COUNTRY date=$DATE short_run=$SHORT_RUN
