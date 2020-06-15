#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
DATE=${1:-$TODAY}

echo "*** Date: $DATE"

echo "*** Copying results"
./run/copy_meffs.R $DATE

