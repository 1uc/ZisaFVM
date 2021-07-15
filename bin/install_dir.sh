#! /usr/bin/env bash

# SPDX-License-Identifier: MIT
# Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

set -e

if [[ "$#" -lt 2 ]]
then
    echo "Usage: $0 COMPILER DESTINATION [--zisa_has_mpi=ZISA_HAS_MPI]"
    echo "                               [--zisa_has_cuda=ZISA_HAS_CUDA]"
    exit -1
fi

for arg in "$@"
do
    case $arg in
        --zisa_has_mpi=*)
            ZISA_HAS_MPI=${arg#*=}
            ;;
        --zisa_has_cuda=*)
            ZISA_HAS_CUDA=${arg#*=}
            ;;
        *)
            ;;
    esac
done

if [[ -z "${ZISA_HAS_MPI}" ]]
then
    ZISA_HAS_MPI=0
fi

if [[ -z "${ZISA_HAS_CUDA}" ]]
then
    ZISA_HAS_CUDA=0
fi

compiler="$1"
compiler_id="$(basename "${compiler}")"
compiler_version="$("$compiler" -dumpversion)"

build_mangled_name="${compiler}__${compiler_version}__ZISA_HAS_MPI=${ZISA_HAS_MPI}__ZISA_HAS_CUDA=${ZISA_HAS_CUDA}"
build_sha="$(echo "${build_mangled_name}" | sha1sum | head -c 6)"

install_root=$(realpath -ms "$2")
install_dir="${install_root}/${build_sha}"

echo "${install_dir}"
