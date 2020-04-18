#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

DATE=$(date "+%Y-%m-%d")

echo "*** ECDC data"
./orderly run ecdc date=$DATE

# Parallel
grep -E '^[A-Z]{3}\s*' countries | \
    parallel -j 4 ./orderly run lmic_reports iso3c={} date=$DATE

# Serial (useful if debugging)
# for ISO in $(grep -E '^[A-Z]{3}\s*' countries); do
#     echo "*** $ISO"
#     ./orderly run lmic_reports iso3c=$ISO date=$DATE
# done

echo "*** Copying reports"
./copy_reports.R $DATE

echo "*** Index page"
./orderly run index_page date=$DATE
echo "*** Parameters page"
./orderly run parameters date=$DATE

echo "*** Copying files"
./copy_index.R $DATE
