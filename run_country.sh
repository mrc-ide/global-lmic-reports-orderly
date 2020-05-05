#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
COUNTRY=$1
DATE=${2:-$TODAY}

echo "*** Country: $COUNTRY"
echo "*** Date: $DATE"

# Parallel
./orderly run lmic_reports iso3c=$COUNTRY date=$DATE
