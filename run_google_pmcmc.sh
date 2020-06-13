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

DEFAULT_HICs="FALSE"
HICs=${5:-$DEFAULT_HICs}

echo "*** Date: $DATE"

echo "*** Short Run: $SHORT_RUN"

echo "*** Full Scenarios: $FULL_SCENARIOS"

echo "*** HICs: $HICs"

echo "*** ECDC data"
./orderly run ecdc date=$DATE

echo "*** Oxford GRT data"
./orderly run oxford_grt date=$DATE

echo "*** Google BRT data"
./orderly run brt_google_mobility date=$DATE short_run=$SHORT_RUN

echo "*** Updating country list"
./update_run_sh.R date=$DATE HICSs=$HICs

echo "*** Running country reports"

# Parallel
grep -E '^[A-Z]{3}\s*' countries | \
parallel --progress -j 64 ./orderly run lmic_reports_google_pmcmc iso3c={} date=$DATE short_run=$SHORT_RUN parallel=$PARALLEL full_scenarios=$FULL_SCENARIOS

# Serial (useful if debugging)
# for ISO in $(grep -E '^[A-Z]{3}\s*' countries); do
#     echo "*** `- $ISO"
#     ./orderly run lmic_reports iso3c=$ISO date=$DATE
# done

echo "*** Copying reports"
./copy_reports_google_pmcmc.R $DATE

echo "*** Index page"
./orderly run index_page date=$DATE

echo "*** Africa page"
./orderly run regional_page date=$DATE continent=Africa
echo "*** Asia page"
./orderly run regional_page date=$DATE continent=Asia
echo "*** Americas page"
./orderly run regional_page date=$DATE continent=Americas
echo "*** Europe page"
./orderly run regional_page date=$DATE continent=Europe

echo "*** Parameters page"
./orderly run parameters date=$DATE
echo "*** 404 page"
./orderly run 404 date=$DATE
echo "*** FAQ page"
./orderly run FAQ date=$DATE
echo "*** News page"
./orderly run news date=$DATE

echo "*** data schema"
./write_data_schema.R

echo "*** Copying files"
./copy_index.R $DATE
./copy_regionals.R $DATE
