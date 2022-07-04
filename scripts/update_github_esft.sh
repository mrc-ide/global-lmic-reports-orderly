#!/usr/bin/env bash
set -e
# ./orderly migrate
# ./orderly rebuild

TODAY=$(date "+%Y-%m-%d")
DATE=${1:-$TODAY}

echo "*** Deploy Key"
REMOTE_URL=git@github.com:mrc-ide/global_lmic_projections_esft.git

export GIT_SSH_COMMAND="ssh -i ../.ssh/gh_esft/id_ed25519"

git -C gh-esft config user.email "gregbarnsley@hotmail.co.uk"
git -C gh-esft config user.name "GBarnsley"

echo "*** Create Commits"
git -C gh-esft add .
git -C gh-esft commit -m "Update projections for version ${DATE}"

echo "*** Push to GitHub"
git -C gh-esft push
