#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
DATE=${1:-$TODAY}

echo "*** Date: $DATE"

echo "*** Copying reports"
./run/copy_reports_google_pmcmc.R $DATE

echo "*** Index page"
./orderly run index_page date=$DATE

echo "*** Africa page"
./orderly run regional_page date=$DATE continent=Africa
echo "*** Asia page"
./orderly run regional_page date=$DATE continent=Asia
echo "*** Americas page"
./orderly run regional_page date=$DATE continent=Americas
echo "*** Europe page"
./orderly run regional_page date=$DATE continent=Europe

echo "*** Parameters page"
./orderly run parameters date=$DATE
echo "*** 404 page"
./orderly run 404 date=$DATE
echo "*** FAQ page"
./orderly run FAQ date=$DATE
echo "*** News page"
./orderly run news date=$DATE

echo "*** data schema"
./run/write_data_schema.R

echo "*** Copying files"
./run/copy_index.R $DATE
./run/copy_regionals.R $DATE