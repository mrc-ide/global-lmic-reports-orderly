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

DEFAULT_HOSP_BEDS="10000000000"
HOSP_BEDS=${4:-$DEFAULT_HOSP_BEDS}

DEFAULT_ICU_BEDS="10000000000"
ICU_BEDS=${5:-$DEFAULT_ICU_BEDS}

echo "*** RF: $RF"
echo "*** Replicates: $REPLICATES"
echo "*** MCMC Iterations: $N_MCMC"
echo "*** Hospital Beds: $HOSP_BEDS"
echo "*** ICU Beds: $ICU_BEDS"

# First batch of 18
parallel -k echo ::: "Andaman and Nicobar Islands" "Andhra Pradesh" "Arunachal Pradesh" "Assam" "Bihar" "Chandigarh" "Chhattisgarh" "Dadra and Nagar Haveli and Daman and Diu" "Delhi" "Goa" "Gujarat" "Haryana" "Himachal Pradesh" "Jammu and Kashmir" "Jharkhand" "Karnataka" "Kerala" "Ladakh"> prov-file

# Parallel
grep -E "*." prov-file | \
parallel --progress -j 18 ./orderly run india_sub_national state={} rf=$RF replicates=$REPLICATES n_mcmc=$N_MCMC hosp_beds=$HOSP_BEDS icu_beds=$ICU_BEDS

rm prov-file

# Second batch of 18
parallel -k echo ::: "Lakshadweep" "Madhya Pradesh" "Maharashtra" "Manipur" "Meghalaya" "Mizoram" "Nagaland" "Odisha" "Puducherry" "Punjab" "Rajasthan" "Sikkim" "Tamil Nadu" "Telangana" "Tripura" "Uttar Pradesh" "Uttarakhand" "West Bengal" > prov-file

# Parallel
grep -E "*." prov-file | \
parallel --progress -j 18 ./orderly run india_sub_national state={} rf=$RF replicates=$REPLICATES n_mcmc=$N_MCMC hosp_beds=$HOSP_BEDS icu_beds=$ICU_BEDS

rm prov-file