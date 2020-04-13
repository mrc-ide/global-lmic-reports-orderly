#!/usr/bin/env bash
set -e
./orderly run ecdc date=$(date "+%Y-%m-%d")

## Fill these in for now manually. Eventually we will use the batch
## interface.
./orderly run lmic_reports country=Angola
./orderly run lmic_reports country=Senegal
