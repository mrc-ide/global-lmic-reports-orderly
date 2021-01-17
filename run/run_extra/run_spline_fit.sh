#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

DEFAULT_RF="1"
rf=${1:-$DEFAULT_RF}

DEFAULT_COUNTRY="Syria"
country=${2:-$DEFAULT_COUNTRY}

DEFAULT_MCMC="100000"
n_mcmc=${3:-$DEFAULT_MCMC}

DEFAULT_HOSP_BEDS="1000000000"
hosp_beds=${4:-$DEFAULT_HOSP_BEDS}

DEFAULT_ICU_BEDS="1000000000"
icu_beds=${5:-$DEFAULT_ICU_BEDS}

echo "*** Reporting Fraction: $rf"
echo "*** Country: $country"
echo "*** MCMC: $n_mcmc"
echo "*** Hospital Beds: $hosp_beds"
echo "*** ICU Beds: $icu_beds"

# Orderly Run
./orderly run spline_fit rf=$rf country=$country n_mcmc=$n_mcmc hosp_beds=$hosp_beds icu_beds=$icu_beds
