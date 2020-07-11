#!/usr/bin/env bash

if (( $# < 1 )) || (( 2 < $# )) ; then
    echo "Usage: $0 COMPILER [LOCATION]"
    exit -1
fi

CC=$1
CC_VER=$(${CC} -dumpversion)

if (( $# == 2 )) ; then
    LOCAL_URL=$(realpath $2)
fi

GIT_URL="gitlab.math.ethz.ch:karoger/helmholts_eos.git"
BUILD_DIR=build-helmholtz
INSTALL_DIR="../../../third_party/helmholtz_eos/${CC}/${CC_VER}"

rm -rf $BUILD_DIR && mkdir $BUILD_DIR && cd $BUILD_DIR

if (( $# == 1 )) ; then
    git clone git@${GIT_URL}
    if (( $? != 0 )); then
        echo "Is VPN running?"
        exit -1
    fi
else
    cp -r ${LOCAL_URL} .
fi

cd helmholts_eos

mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} ..
make
make install
