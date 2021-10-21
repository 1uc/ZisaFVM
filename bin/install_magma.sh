#! /usr/bin/env bash

set -e
set -o xtrace

if [[ $# -ne 1 ]]; then
  echo "Usage:"
  echo "    ${0} INSTALL_PREFIX"
  exit -1
fi

INSTALL_PREFIX=$(realpath $1)

MAGMA_VERSION=2.6.1



TMP_DIR=$(mktemp -d --tmpdir=${SCRATCH})
MAGMA_DIR="magma-${MAGMA_VERSION}"
URL="http://icl.utk.edu/projectsfiles/magma/downloads/magma-${MAGMA_VERSION}.tar.gz"

pushd $TMP_DIR
wget ${URL}
tar xf magma-${MAGMA_VERSION}.tar.gz
pushd ${MAGMA_DIR}

read -r -d '' make_inc_patch <<'EOF' || true
--- ../../tmp.cztUVvsu15/magma-2.6.1/make.inc	2021-07-13 00:35:20.000000000 +0200
+++ make.inc	2021-10-16 15:27:30.650000209 +0200
@@ -55,8 +55,7 @@
 # set our GPU targets
 ifeq ($(BACKEND),cuda)
     # For newer cards
-    #GPU_TARGET = Pascal Volta Turing Ampere
-    GPU_TARGET = Kepler Maxwell Pascal
+    GPU_TARGET = Maxwell Pascal Volta Turing Ampere
 else ifeq ($(BACKEND),hip)
     GPU_TARGET = gfx803 gfx900 gfx901
 endif
EOF
cp make.inc-examples/make.inc.openblas make.inc
patch make.inc <(echo "${make_inc_patch}")

export OPENBLASDIR=/usr
export CUDADIR=/opt/cuda

make -j$(nproc)
make install prefix=${INSTALL_PREFIX}/${MAGMA_DIR}

# cd ${HOME}
# rm -rf $TMP_DIR
