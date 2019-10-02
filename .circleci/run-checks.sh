#!/bin/bash

set -e
set -eo pipefail

source activate ${ENV_NAME}
${NM:-nm} ${CONDA_PREFIX}/lib/libncomp.a
if [ $(uname) = Darwin ]; then
    ${OTOOL:-otool} -L ${CONDA_PREFIX}/lib/libncomp.dylib
else
    ${LDD:-ldd} ${CONDA_PREFIX}/lib/libncomp.so
fi
make check
