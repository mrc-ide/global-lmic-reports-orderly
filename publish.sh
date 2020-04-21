#!/usr/bin/env bash
set -e
DOCS_DIR=gh-pages
VERSION=$(git rev-parse --short HEAD)

WHERE=$1

if [[ $WHERE == "staging" ]]; then
    REMOTE_URL=git@github.com:mrc-ide/global-lmic-reports-staging.git
elif [[ $WHERE == "production" ]]; then
    REMOTE_URL=git@github.com:mrc-ide/global-lmic-reports.git
else
    echo "Usage: ./publish.sh staging|production"
    exit 1
fi

if [ -f .ssh/id_rsa ]; then
    echo "Using deploy key"
    ## NOTE, uses one directory above the root
    export GIT_SSH_COMMAND="ssh -i ../.ssh/id_rsa"
fi

rm -rf ${DOCS_DIR}/.git
git init ${DOCS_DIR}
git -C ${DOCS_DIR} add .
git -C ${DOCS_DIR} commit --no-verify -m "Update pages for version ${VERSION}"
git -C ${DOCS_DIR} remote add origin ${REMOTE_URL}
git -C ${DOCS_DIR} push --force -u origin master
