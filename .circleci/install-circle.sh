#!/bin/bash

set -e
set -eo pipefail

conda config --set always_yes true --set changeps1 false --set quiet true
conda config --add channels conda-forge
conda env create -f .circleci/environment-dev-$(uname).yml --name=${ENV_NAME} --quiet
conda env list
source activate ${ENV_NAME}
conda install -y -c conda-forge blas
conda install -y -c conda-forge lapack
autoreconf --install
./configure --prefix=${CONDA_PREFIX}
make install
