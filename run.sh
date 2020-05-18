#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
DATE=${1:-$TODAY}
DEFAULT_SHORT="FALSE"
SHORT_RUN=${2:-$DEFAULT_SHORT}

echo "*** Date: $DATE"

echo "*** ECDC data"
./orderly run ecdc date=$DATE

echo "*** Oxford GRT data"
./orderly run oxford_grt date=$DATE

# echo "*** Google BRT data"
# ./orderly run brt_google_mobility date=$DATE

echo "*** Updating country list"
./update_run_sh.R $DATE

echo "*** Running country reports"

# Parallel
grep -E '^[A-Z]{3}\s*' countries | \
    parallel --progress -j 32 ./orderly run lmic_reports iso3c={} date=$DATE short_run=$SHORT_RUN
# Serial (useful if debugging)
# for ISO in $(grep -E '^[A-Z]{3}\s*' countries); do
#     echo "*** `- $ISO"
#     ./orderly run lmic_reports iso3c=$ISO date=$DATE
# done

echo "*** Copying reports"
./copy_reports.R $DATE

echo "*** Index page"
./orderly run index_page date=$DATE
echo "*** Parameters page"
./orderly run parameters date=$DATE
echo "*** 404 page"
./orderly run 404 date=$DATE
echo "*** FAQ page"
./orderly run FAQ date=$DATE

echo "*** data schema"
./write_data_schema.R

echo "*** Copying files"
./copy_index.R $DATE
