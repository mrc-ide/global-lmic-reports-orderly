#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

DEFAULT_RF="1"
RF=${1:-$DEFAULT_RF}

DEFAULT_REPLICATES="20"
REPLICATES=${2:-$DEFAULT_REPLICATES}

DEFAULT_N_MCMC="100000"
N_MCMC=${3:-$DEFAULT_N_MCMC}

DEFAULT_PROB_HOSP_MULTIPLIER="1"
PROB_HOSP_MULTIPLIER=${4:-$DEFAULT_PROB_HOSP_MULTIPLIER}

DEFAULT_DUR_R="365"
DUR_R=${5:-$DEFAULT_DUR_R}

DEFAULT_MODEL="SQUIRE"
MODEL=${6:-$DEFAULT_MODEL}

DEFAULT_STATE="Kerala"
STATE=${7:-$DEFAULT_STATE}

TODAY=$(date "+%Y-%m-%d")
DATE=${8:-$TODAY}


echo "*** RF: $RF"
echo "*** Replicates: $REPLICATES"
echo "*** MCMC Iterations: $N_MCMC"
echo "*** Prob Hosp Multiplier: $PROB_HOSP_MULTIPLIER"
echo "*** Dur R: $DUR_R"
echo "*** Model: $MODEL"
echo "*** State: $STATE"

# Parallel
parallel --progress -j 18 ./orderly run india_sub_national state=$STATE rf=$RF replicates=$REPLICATES n_mcmc=$N_MCMC prob_hosp_multiplier=$PROB_HOSP_MULTIPLIER dur_R=$DUR_R model=$MODEL date=$DATE
