#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild


COUNTRY=$1

TODAY=$(date "+%Y-%m-%d")
DATE=${2:-$TODAY}

DEFAULT_SHORT="FALSE"
SHORT_RUN=${3:-$DEFAULT_SHORT}

DEFAULT_PARALLEL="TRUE"
PARALLEL=${4:-$DEFAULT_PARALLEL}

DEFAULT_FULL_SCENARIOS="FALSE"
FULL_SCENARIOS=${5:-$DEFAULT_FULL_SCENARIOS}

DEFAULT_GIBBS="FALSE"
GIBBS=${6:-$DEFAULT_GIBBS}

DEFAULT_DUR_R=365
DUR_R=${7:-$DEFAULT_DUR_R}

DEFAULT_DUR_V=1095
DUR_V=${8:-$DEFAULT_DUR_V}

DEFAULT_VACCINE_UPTAKE=0.95
VACCINE_UPTAKE=${9:-$DEFAULT_VACCINE_UPTAKE}

DEFAULT_AVD=0.95
AVD=${10:-$DEFAULT_AVD}

DEFAULT_VACINE_SCENARIO="UK"
VACCINE_SCENARIO=${11:-$DEFAULT_VACINE_SCENARIO}


echo "*** Country: $COUNTRY"
echo "*** Date: $DATE"
echo "*** Short Run: $SHORT_RUN"
echo "*** Parallel: $PARALLEL"
echo "*** Full Scenarios: $FULL_SCENARIOS"
echo "*** GIBBS: $GIBBS"

echo "*** dur_R: $DUR_R"
echo "*** dur_V: $DUR_V"
echo "*** Vaccine Uptake: $VACCINE_UPTAKE"
echo "*** Available Doses Proportion: $DEFAULT_AVD"
echo "*** Vaccine Scenario: $VACCINE_SCENARIO"


# Run
./orderly run vaccine_control_fit iso3c=$COUNTRY \
date=$DATE short_run=$SHORT_RUN parallel=$PARALLEL \
full_scenarios=$FULL_SCENARIOS gibbs_sampling=$GIBBS \
dur_R=$DUR_R dur_V=$DUR_V vaccine_uptake=$VACCINE_UPTAKE \
available_doses_proportion=$AVD vaccine_scenario=$VACCINE_SCENARIO