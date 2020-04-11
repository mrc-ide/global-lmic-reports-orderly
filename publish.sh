#!/bin/sh
set -e
DOCS_DIR=gh-pages
VERSION=$(git rev-parse --short HEAD)
REMOTE_URL=git@github.com:mrc-ide/global-lmic-reports.git

if [ -d $DOCS_DIR ]; then
    git -C $DOCS_DIR pull
else
    git clone $REMOTE_URL $DOCS_DIR
fi

rm -rf ${DOCS_DIR}/.git
git init ${DOCS_DIR}
git -C ${DOCS_DIR} add .
git -C ${DOCS_DIR} commit --no-verify -m "Update pages for version ${VERSION}"
git -C ${DOCS_DIR} remote add origin ${REMOTE_URL}
git -C ${DOCS_DIR} push --force -u origin master
