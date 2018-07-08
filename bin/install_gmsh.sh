#! /usr/bin/env bash
set -x

TOP="${PWD}"

GMSH_VERSION="3.0.6"
GMSH_DIR=third_party/gmsh-${GMSH_VERSION}
GMSH_BIN_DIR=${GMSH_DIR}/bin
GMSH_NAME=gmsh-${GMSH_VERSION}-source
GMSH_TAR_NAME=${GMSH_NAME}.tgz

mkdir -p ${GMSH_BIN_DIR}

wget http://gmsh.info/src/gmsh-${GMSH_VERSION}-source.tgz -O ${GMSH_DIR}/${GMSH_TAR_NAME}

cd ${GMSH_DIR}
echo ${PWD}
tar xf ${GMSH_TAR_NAME}

cd ${GMSH_NAME}
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../.. -DCMAKE_BUILD_TYPE=Release ..
make -j$(($(nproc) + 1)) && make install

ln -sf ${TOP}/${GMSH_BIN_DIR}/gmsh ${TOP}/bin/gmsh
