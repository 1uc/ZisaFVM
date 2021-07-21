#!/usr/bin/env bash

# SPDX-License-Identifier: MIT
# Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

set -e

if (( $# != 1 )); then
  echo "Usage: $0 [compiler]"
  exit -1
fi

CC=$1
CC_VER=$(${CC} -dumpversion)

BUILD_DIR="build-third_party"

METIS_VER="5.1.0"
METIS_LINK="http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-${METIS_VER}.tar.gz"
METIS_DOC_LINK="http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf"

METIS_INSTALL="${PWD}/third_party/metis-${METIS_VER}/${CC}/${CC_VER}"

mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

if [[ -f metis-${METIS_VER} ]]; then
    rm -rf metis-${METIS_VER}
fi

if [[ ! -f metis-${METIS_VER}.tar.gz ]]; then
  wget ${METIS_LINK}
fi
tar xf metis-${METIS_VER}.tar.gz

cd metis-${METIS_VER}

make config cc=${CC} prefix=${METIS_INSTALL}
make
make install

mkdir -p ${METIS_INSTALL}/doc
wget ${METIS_DOC_LINK} -O ${METIS_INSTALL}/doc/manual.pdf
