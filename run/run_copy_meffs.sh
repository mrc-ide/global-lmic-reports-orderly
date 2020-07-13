#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
DATE=${1:-$TODAY}
DEFAULT_WHICH="BOTH"
WHAT=${2:-$DEFAULT_WHAT}
DEFAULT_DIC_ONLY="FALSE"
DIC_ONLY=${3:-$DEFAULT_DIC}

echo "*** Date: $DATE"
echo "*** Which: $WHAT"
echo "*** DIC Only: $DIC_ONLY"

echo "*** Copying results"
./run/copy_meffs.R $DATE $WHAT $DIC_ONLY

