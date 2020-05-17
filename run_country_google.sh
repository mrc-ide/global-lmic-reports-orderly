#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
COUNTRY=$1
DATE=${2:-$TODAY}
DEFAULT_SHORT="FALSE"
SHORT_RUN=${3:-$DEFAULT_SHORT}
DEFAULT_PARALLEL="FALSE"
PARALLEL=${4:-$DEFAULT_PARALLEL}

echo "*** Country: $COUNTRY"
echo "*** Date: $DATE"
echo "*** Short Run: $SHORT_RUN"
echo "*** Parallel: $PARALLEL"

# Parallel
./orderly run lmic_reports_google iso3c=$COUNTRY date=$DATE short_run=$SHORT_RUN parallel=$PARALLEL