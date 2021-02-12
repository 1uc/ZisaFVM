#! /usr/bin/env bash
set -e

component_name="ZisaFVM"

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
    ZISA_HAS_CUDA=1
fi

if [[ ${ZISA_HAS_MPI} -eq 0 ]]
then
    zisa_dependencies=("ZisaCore" "ZisaMemory" "ZisaSFC" "ZisaTimeStepping")
else
    zisa_dependencies=("ZisaCore" "ZisaMemory" "ZisaSFC" "ZisaMPI" "ZisaTimeStepping")
fi


zisa_root="$(realpath "$(dirname "$(readlink -f "$0")")"/..)"

CC="$1"
CXX="$(${zisa_root}/bin/cc2cxx.sh $CC)"

install_dir="$("${zisa_root}/bin/install_dir.sh" "$1" "$2" --zisa_has_mpi=${ZISA_HAS_MPI})"
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
        cd "${src_dir}"

        if [[ -z $(git remote -v 2>/dev/null | grep ${repo_url}) ]]
        then
            echo "Failed to install ${dep} to ${src_dir}"
            exit -1

        else
            cd "${HOME}"
            rm -rf "${src_dir}"
        fi
    fi

    git clone ${repo_url} "${src_dir}"

    mkdir -p "${src_dir}/build-dep"
    cd "${src_dir}/build-dep"

    cmake -DCMAKE_INSTALL_PREFIX="${install_dir}/zisa" \
          -DCMAKE_PREFIX_PATH="${install_dir}/zisa/lib/cmake/zisa" \
          -DCMAKE_MODULE_PATH="${install_dir}/conan" \
          -DCMAKE_PROJECT_${dep}_INCLUDE="${install_dir}/conan/conan_paths.cmake" \
          -DCMAKE_C_COMPILER="${CC}" \
          -DCMAKE_CXX_COMPILER="${CXX}" \
          -DZISA_HAS_MPI="${ZISA_HAS_MPI}" \
          -DCMAKE_BUILD_TYPE="Release" \
          ..

    cmake --build . --parallel $(nproc)
    cmake --install .
done

if [[ ${ZISA_HAS_METIS} -ne 0 ]]
then
    cd "${zisa_fvm_root}"
    bin/install_metis.sh "$1"
    cd -
fi


echo "The dependencies were installed at"
echo "    export DEP_DIR=${install_dir}"
echo ""
echo "Use"
echo "    cmake -DCMAKE_PROJECT_${component_name}_INCLUDE=${install_dir}/conan/conan_paths.cmake \\ "
echo "          -DCMAKE_MODULE_PATH=${install_dir}/conan \\ "
echo "          -DCMAKE_PREFIX_PATH=${install_dir}/zisa/lib/cmake/zisa \\ "
echo "          -DCMAKE_C_COMPILER=${CC} \\ "
echo "          -DCMAKE_CXX_COMPILER=${CXX} \\ "
echo "          -DZISA_HAS_MPI=${ZISA_HAS_MPI} \\ "
echo "          REMAINING_ARGS "
