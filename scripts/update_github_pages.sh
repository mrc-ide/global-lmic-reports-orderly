#!/usr/bin/env bash
set -e
DOCS_DIR=gh-pages
VERSION=$(git rev-parse --short HEAD)

REMOTE_URL=git@github.com:mrc-ide/global-lmic-reports.git

echo "Using deploy key"
## NOTE, uses one directory above the root
export GIT_SSH_COMMAND="ssh -i ../.ssh/gh_pages/id_ed25519"

if [ -z $(git -C ${DOCS_DIR} config --get user.email) ]; then
    git -C ${DOCS_DIR} config user.email "gregbarnsley@hotmail.co.uk"
    git -C ${DOCS_DIR} config user.name "GBarnsley"
fi

git -C ${DOCS_DIR} add data/
git -C ${DOCS_DIR} commit --no-verify -m "Add data for version ${VERSION}"
git -C ${DOCS_DIR} push --force -u origin main
git -C ${DOCS_DIR} add .
git -C ${DOCS_DIR} commit --no-verify -m "Update pages for version ${VERSION}"
git -C ${DOCS_DIR} push --force -u origin main
git -C ${DOCS_DIR} branch master -d
git -C ${DOCS_DIR} branch master
git -C ${DOCS_DIR} checkout master
git -C ${DOCS_DIR} push --force -u origin master
git -C ${DOCS_DIR} checkout main
