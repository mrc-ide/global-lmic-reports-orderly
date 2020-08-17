#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
DATE=${1:-$TODAY}

DEFAULT_URBAN="FALSE"
URBAN=${2:-$DEFAULT_URBAN}

DEFAULT_PHO="TRUE"
PHO=${3:-$DEFAULT_PHO}

DEFAULT_YC="TRUE"
YC=${4:-$DEFAULT_YC}


echo "*** Date: $DATE"
echo "*** Urban: $URBAN"
echo "*** Poorer Health Outcomes: $PHO"
echo "*** Younger Cities: $YC"

# Parallel
grep -E "*." rf.txt | \
parallel --progress --verbose -j 8 ./orderly run syria_under_reporting reporting_fraction={} date=$DATE urban=$URBAN poorer_health_outcomes=$PHO younger_cities=$YC
