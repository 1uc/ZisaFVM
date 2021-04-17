#! /usr/bin/env bash
set -e

if [[ "$#" -ne 1 ]]
then
    echo "Usage: $0 DESTINATION"
    exit -1
fi

TOP=${PWD}

mkdir -p ${TOP}/tmp

GMSH_VERSION="4.7.1"
GMSH_INSTALL_PREFIX="$1/gmsh-${GMSH_VERSION}"
GMSH_DIR=${TOP}/$(mktemp -d --tmpdir=tmp -t gmsh.XXXXXX)
GMSH_NAME=gmsh-${GMSH_VERSION}-source
GMSH_TAR_NAME=${GMSH_NAME}.tgz

wget http://gmsh.info/src/gmsh-${GMSH_VERSION}-source.tgz -O ${GMSH_DIR}/${GMSH_TAR_NAME}

cd ${GMSH_DIR}
tar xf ${GMSH_TAR_NAME}

cd ${GMSH_NAME}
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX="${GMSH_INSTALL_PREFIX}" -DCMAKE_BUILD_TYPE=Release ..
make -j$(($(nproc) + 1)) && make install

if ! command -v gmsh &> /dev/null
then
  echo "WARN: Could not find 'gmsh'. Please ensure that"
  echo "    ${GMSH_INSTALL_PREFIX}/bin"
  echo "is in the PATH."
fi

cd ${TOP}
# rm -rf ${GMSH_DIR}
