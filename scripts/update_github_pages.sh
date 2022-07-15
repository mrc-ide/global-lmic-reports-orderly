#!/usr/bin/env bash
set -e
VERSION=$(git rev-parse --short HEAD)

REMOTE_URL=git@github.com:mrc-ide/global-lmic-reports.git

echo "Using deploy key"
## NOTE, uses one directory above the root
export GIT_SSH_COMMAND="ssh -i ../.ssh/gh_pages/id_ed25519"

git -C gh-pages config user.email "gregbarnsley@hotmail.co.uk"
git -C gh-pages config user.name "GBarnsley"

git -C gh-pages add data/
git -C gh-pages commit --no-verify -m "Add data for version ${VERSION}"
git -C gh-pages push --force -u origin main
git -C gh-pages add .
git -C gh-pages commit --no-verify -m "Update pages for version ${VERSION}"
git -C gh-pages push --force -u origin main
git -C gh-pages branch master -d
git -C gh-pages branch master
git -C gh-pages checkout master
git -C gh-pages push --force -u origin master
git -C gh-pages checkout main
