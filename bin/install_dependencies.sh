#! /usr/bin/env bash

# SPDX-License-Identifier: MIT
# Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

set -e

component_name="ZisaFVM"

if [[ "$#" -lt 2 ]]
then
    echo "Usage: $0 COMPILER DESTINATION [--zisa_has_mpi=ZISA_HAS_MPI]"
    echo "                               [--zisa_has_netcdf=ZISA_HAS_NETCDF]"
    echo "                               [--cmake=CUSTOM_CMAKE_BINARY]"
    echo "                               [--print_install_dir]"
    exit -1
fi

for arg in "$@"
do
    case $arg in
        --zisa_has_mpi=*)
            ZISA_HAS_MPI=${arg#*=}
            ;;
        --zisa_has_netcdf=*)
            ZISA_HAS_NETCDF=${arg#*=}
            ;;
        --cmake=*)
            CMAKE="$(realpath "${arg#*=}")"
            ;;
        --print_install_dir)
            PRINT_INSTALL_PATH=1
            ;;
        *)
            ;;
    esac
done

if [[ -z "${CMAKE}" ]]
then
    CMAKE=cmake
fi

if [[ -z "${ZISA_HAS_MPI}" ]]
then
    ZISA_HAS_MPI=0
fi

if [[ -z "${ZISA_HAS_HDF5}" ]]
then
    # This should be on since ZisaFVM uses HDF5 for I/O.
    ZISA_HAS_HDF5=1
fi

if [[ -z "${ZISA_HAS_NETCDF}" ]]
then
    ZISA_HAS_NETCDF=0
fi

if [[ ${ZISA_HAS_MPI} -eq 0 ]]
then
    zisa_dependencies=("ZisaCore" "ZisaMemory" "ZisaSFC" "ZisaTimeStepping")
else
    zisa_dependencies=("ZisaCore" "ZisaMemory" "ZisaSFC" "ZisaMPI" "ZisaTimeStepping")
fi

if [[ -z "${ZISA_HAS_CUDA}" ]]
then
    # This should be off since ZisaFVM doesn't support CUDA.
    ZISA_HAS_CUDA=0
fi


zisa_root="$(realpath "$(dirname "$(readlink -f "$0")")"/..)"

CC="$1"
CXX="$(${zisa_root}/bin/cc2cxx.sh $CC)"

install_dir="$(
    "${zisa_root}/bin/install_dir.sh" "$1" "$2" \
        --zisa_has_mpi=${ZISA_HAS_MPI} \
        --zisa_has_cuda=${ZISA_HAS_CUDA} \
        --zisa_has_hdf5=${ZISA_HAS_HDF5} \
        --zisa_has_netcdf=${ZISA_HAS_NETCDF} \
)"

if [[ ${PRINT_INSTALL_PATH} -eq 1 ]]
then
  echo $install_dir
  exit 0
fi

source_dir="${install_dir}/sources"
conan_file="${zisa_root}/conanfile.txt"

mkdir -p "${install_dir}/conan" && cd "${install_dir}/conan"
conan install "$conan_file" \
        -s compiler=$(basename "${CC}") \
        -s compiler.libcxx=libstdc++11

mkdir -p "${source_dir}"
for dep in "${zisa_dependencies[@]}"
do
    src_dir="${source_dir}/$dep"
    repo_url=git@github.com:1uc/${dep}.git

    # If necessary and reasonable remove ${src_dir}.
    if [[ -d "${src_dir}" ]]
    then
        pushd "${src_dir}"

        if [[ -z $(git remote -v 2>/dev/null | grep ${repo_url}) ]]
        then
            echo "Failed to install ${dep} to ${src_dir}"
            exit -1

        else
            popd
            rm -rf "${src_dir}"
        fi
    fi

    git clone ${repo_url} "${src_dir}"
    pushd "${src_dir}"

    build_dir="${src_dir}/build"
    ${CMAKE} -DCMAKE_INSTALL_PREFIX="${install_dir}/zisa" \
             -DCMAKE_PREFIX_PATH="${install_dir}/zisa/lib/cmake/zisa" \
             -DCMAKE_PROJECT_${dep}_INCLUDE="${install_dir}/conan/conan_paths.cmake" \
             -DCMAKE_C_COMPILER="${CC}" \
             -DCMAKE_CXX_COMPILER="${CXX}" \
             -DZISA_HAS_CUDA="${ZISA_HAS_CUDA}" \
             -DZISA_HAS_MPI="${ZISA_HAS_MPI}" \
             -DZISA_HAS_HDF5="${ZISA_HAS_HDF5}" \
             -DZISA_HAS_NETCDF="${ZISA_HAS_NETCDF}" \
             -DCMAKE_BUILD_TYPE="Release" \
             -B ${build_dir}

    ${CMAKE} --build ${build_dir} --parallel $(nproc)
    ${CMAKE} --install ${build_dir}

    popd
done

if [[ ${ZISA_HAS_METIS} -ne 0 ]]
then
    cd "${zisa_root}"
    bin/install_metis.sh "$1"
    cd -
fi


echo "The dependencies were installed at"
echo "    export DEP_DIR=${install_dir}"
echo ""
echo "Use"
echo "    ${CMAKE} -DCMAKE_PROJECT_${component_name}_INCLUDE=${install_dir}/conan/conan_paths.cmake \ "
echo "             -DCMAKE_PREFIX_PATH=${install_dir}/zisa/lib/cmake/zisa \ "
echo "             -DCMAKE_C_COMPILER=${CC} \ "
echo "             -DCMAKE_CXX_COMPILER=${CXX} \ "
echo "             -DCMAKE_BUILD_TYPE=FastDebug \ "
echo "             -B build"
