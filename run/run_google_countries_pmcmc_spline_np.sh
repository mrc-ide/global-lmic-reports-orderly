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

DEFAULT_PARALLEL="TRUE"
PARALLEL=${4:-$DEFAULT_PARALLEL}

DEFAULT_HICs="TRUE"
HICs=${5:-$DEFAULT_HICs}

DEFAULT_COUNTRIES="countries"
COUNTRIES=${6:-$DEFAULT_COUNTRIES}

echo "*** Date: $DATE"

echo "*** Short Run: $SHORT_RUN"

echo "*** Full Scenarios: $FULL_SCENARIOS"

echo "*** Parallel: $PARALLEL"

echo "*** HICs: $HICs"

echo "*** Countries: $COUNTRIES"

#echo "*** Updating country list"
#./run/update_run_sh.R $DATE $HICs

echo "*** Running country reports"

# Parallel
grep -E '^[A-Z]{3}\s*' $COUNTRIES | \
parallel --progress -j 62  ./orderly run lmic_reports_google_pmcmc_spline_np iso3c={} date=$DATE short_run=$SHORT_RUN parallel=$PARALLEL full_scenarios=$FULL_SCENARIOS

