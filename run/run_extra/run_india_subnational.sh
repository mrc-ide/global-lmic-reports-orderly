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

TODAY=$(date "+%Y-%m-%d")
DATE=${7:-$TODAY}


echo "*** RF: $RF"
echo "*** Replicates: $REPLICATES"
echo "*** MCMC Iterations: $N_MCMC"
echo "*** Prob Hosp Multiplier: $PROB_HOSP_MULTIPLIER"
echo "*** Dur R: $DUR_R"
echo "*** Model: $MODEL"

# Batch in 18
parallel -k echo ::: "Bihar" "Uttar Pradesh" "Tamil Nadu" "Maharashtra" "Odisha" "West Bengal" "Telangana" "Karnataka" "Andhra Pradesh" "Madhya Pradesh" "Assam" "Rajasthan" "Jharkhand" "Chhattisgarh" "Delhi" "Punjab" "Jammu and Kashmir" "Kerala" "Gujarat" "Haryana" "Uttarakhand" "Himachal Pradesh" "Tripura" "Manipur" "Nagaland" "Goa" "Puducherry" "Meghalaya" "Chandigarh" "Arunachal Pradesh" "Andaman and Nicobar Islands" "Sikkim" "Ladakh" "Mizoram" "Lakshadweep" "Dadra and Nagar Haveli and Daman and Diu"> prov-file

# Parallel
grep -E "*." prov-file | \
parallel --progress -j 18 ./orderly run india_sub_national state={} rf=$RF replicates=$REPLICATES n_mcmc=$N_MCMC prob_hosp_multiplier=$PROB_HOSP_MULTIPLIER dur_R=$DUR_R model=$MODEL date=$DATE

rm prov-file
