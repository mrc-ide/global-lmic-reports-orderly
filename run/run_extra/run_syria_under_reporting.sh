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

DEFAULT_CA="older"
CA=${4:-$DEFAULT_CA}

DEFAULT_HNU="0.2"
HNU=${5:-$DEFAULT_HNU}

echo "*** Date: $DATE"
echo "*** Urban: $URBAN"
echo "*** Poorer Health Outcomes: $PHO"
echo "*** City Age: $CA"
echo "*** Hospital Normal Use: $HNU"

# Parallel
grep -E "*." rf.txt | \
parallel --progress -j 8 ./orderly run syria_under_reporting reporting_fraction={} date=$DATE urban=$URBAN poorer_health_outcomes=$PHO city_age=$CA hospital_normal_use=$HNU
