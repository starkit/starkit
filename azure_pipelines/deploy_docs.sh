#!/bin/bash

# Build the documentation from the SOURCE_BRANCH
# and push it to TARGET_BRANCH.
SOURCE_BRANCH="master"
TARGET_BRANCH="gh-pages"


# Store some useful information
REPO=`git config remote.origin.url`
SSH_REPO=${REPO/https:\/\/github.com\//git@github.com:}
SHA=`git rev-parse --verify HEAD` # Latest SHA Commit

# SOURCE_BRANCH (master) can be acessed from default dir at VM
# To access TARGET_BRANCH (gh-pages) alongside it, we use out/ dir
git clone $REPO out
cd out
# If gh-pages doesn't exist yet (should only happen on first deploy),
# create a new empty (orphan) branch 
git checkout $TARGET_BRANCH || git checkout --orphan $TARGET_BRANCH
# Clean out existing contents
git rm -rf . || exit 0
# Move back to SOURCE_BRANCH
cd ..


# Build the Sphinx documentation
python setup.py build_docs
# Move built docs to out/
mv -f docs/_build/html/* out/
touch out/.nojekyll


# Move again to out/ (TARGET_BRANCH)
cd out
git config --local user.name "Azure Pipelines"
git config --local user.email "azuredevops@microsoft.com"

echo "Doing git add/commit/push"
# Stage all the "changes" i.e. the new version
git add --all

# Exit if there are no difference between new & older version
if git diff --staged --quiet; then
  echo "Exiting with no docs changes"
  exit 0
fi

# Otherwise, commit and push
git commit -m "Deploy to GitHub Pages: ${SHA}"
git push $SSH_REPO $TARGET_BRANCH
cd ..
