#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
DATE=${1:-$TODAY}

DEFAULT_POP="2079000"
DAM_POP=${2:-$DEFAULT_POP}

DEFAULT_PHO="TRUE"
PHO=${3:-$DEFAULT_PHO}

DEFAULT_CA="older"
CA=${4:-$DEFAULT_CA}

DEFAULT_HOSP="2000"
HOSP=${5:-$DEFAULT_HOSP}

DEFAULT_HNU="0.2"
HNU=${6:-$DEFAULT_HNU}

DEFAULT_LS="FALSE"
LS=${7:-$DEFAULT_LS}

DEFAULT_DTF="reported"
DTF=${8:-$DEFAULT_DTF}

echo "*** Date: $DATE"
echo "*** Urban: $URBAN"
echo "*** Poorer Health Outcomes: $PHO"
echo "*** City Age: $CA"
echo "*** Hospital Beds: $HOSP"
echo "*** Hospital Normal Use: $HNU"
echo "*** Late Start: $LS"
echo "*** Data To Fit: $DTF"

# Parallel
grep -E "*." rf.txt | \
parallel --progress -j 11 ./orderly run syria_under_reporting reporting_fraction={} date=$DATE dam_pop=$DAM_POP poorer_health_outcomes=$PHO city_age=$CA hosp_beds=$HOSP hospital_normal_use=$HNU late_start=$LS data_to_fit=$DTF 
