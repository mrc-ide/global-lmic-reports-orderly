#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
DATE=${1:-$TODAY}
DEFAULT_SHORT="FALSE"
SHORT_RUN=${2:-$DEFAULT_SHORT}
DEFAULT_FULL_SCENARIOS="FALSE"
FULL_SCENARIOS=${3:-$DEFAULT_FULL_SCENARIOS}

echo "*** Date: $DATE"

echo "*** Short Run: $SHORT_RUN"

echo "*** Full Scenarios: $FULL_SCENARIOS"

echo "*** ECDC data"
./orderly run ecdc date=$DATE

echo "*** Oxford GRT data"
./orderly run oxford_grt date=$DATE

echo "*** Google BRT data"
./orderly run brt_google_mobility date=$DATE short_run=$SHORT_RUN