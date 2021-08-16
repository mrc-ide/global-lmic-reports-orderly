#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
COUNTRY=$1
DATE=${2:-$TODAY}
DEFAULT_SHORT="FALSE"
SHORT_RUN=${3:-$DEFAULT_SHORT}
DEFAULT_PARALLEL="TRUE"
PARALLEL=${4:-$DEFAULT_PARALLEL}
DEFAULT_FULL_SCENARIOS="FALSE"
FULL_SCENARIOS=${5:-$DEFAULT_FULL_SCENARIOS}
DEFAULT_GIBBS="FALSE"
GIBBS=${6:-$DEFAULT_GIBBS}
DEFAULT_N_MCMC=20000
N_MCMC=${6:-$DEFAULT_N_MCMC}

echo "*** Country: $COUNTRY"
echo "*** Date: $DATE"
echo "*** Short Run: $SHORT_RUN"
echo "*** Parallel: $PARALLEL"
echo "*** Full Scenarios: $FULL_SCENARIOS"
echo "*** GIBBS: $GIBBS"
echo "*** N_MCMC: $N_MCMC"

# Parallel
./orderly run lmic_reports_vaccine iso3c=$COUNTRY \
date=$DATE short_run=$SHORT_RUN parallel=$PARALLEL \
full_scenarios=$FULL_SCENARIOS gibbs_sampling=$GIBBS \
n_mcmc=$N_MCMC
