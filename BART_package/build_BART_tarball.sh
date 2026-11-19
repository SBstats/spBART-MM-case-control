#!/bin/bash
################################################################################
# build_BART_tarball.sh
#
# Rebuild BART_2.9.9.1.tar.gz from CRAN BART 2.9.9 plus the files in
# modified_files/. Only needed to rebuild the tarball after changing modified_files/.
#
# Usage (from anywhere):  bash BART_package/build_BART_tarball.sh
################################################################################
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

cd "$WORK"
curl -sfL -o BART_2.9.9.tar.gz https://cran.r-project.org/src/contrib/Archive/BART/BART_2.9.9.tar.gz \
  || curl -sfL -o BART_2.9.9.tar.gz https://cran.r-project.org/src/contrib/BART_2.9.9.tar.gz
tar -xzf BART_2.9.9.tar.gz
rm -f BART/MD5

cp "$HERE/modified_files/DESCRIPTION" "$HERE/modified_files/NAMESPACE" BART/
cp "$HERE/modified_files/R/wbart_stateful.R" BART/R/
cp "$HERE/modified_files/src/cwbart_stateful.cpp" "$HERE/modified_files/src/init.c" BART/src/

R CMD build --no-build-vignettes --no-manual BART
VERSION="$(grep '^Version:' BART/DESCRIPTION | awk '{print $2}')"
cp "BART_${VERSION}.tar.gz" "$HERE/"
echo "Built $HERE/BART_${VERSION}.tar.gz"
echo "If the version changed, update SPBART_BART_VERSION in R/spBART_sampler.R."
