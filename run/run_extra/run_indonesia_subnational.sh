#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

DEFAULT_ALL_DEATHS="TRUE"
ALL_DEATHS=${1:-$DEFAULT_ALL_DEATHS}

echo "*** Province: $PROVINCE"
echo "*** All Deaths: $ALL_DEATHS"

parallel -k echo ::: "BANTEN" "CENTRAL JAVA" "EAST JAVA"  "JAKARTA"  "WEST JAVA" "YOGYAKARTA" > prov-file

# Parallel
grep -E "*." prov-file | \
parallel --progress -j 11 ./orderly run indonesia_sub_national province={} all_deaths=$ALL_DEATHS

rm prov-file