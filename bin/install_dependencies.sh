#! /usr/bin/env bash

if [[ "$#" -lt 2 ]]
then
    echo "Usage: $0 COMPILER DESTINATION [--zisa_has_mpi=ZISA_HAS_MPI]"
fi

for arg in "$@"
do
    case $arg in
        --zisa_has_mpi=*)
            ZISA_HAS_MPI=${arg#*=}
            ;;
        *)
            ;;
    esac
done

if [[ -z "${ZISA_HAS_MPI}" ]]
then
    ZISA_HAS_MPI=0
fi

component_name="ZisaFVM"
if [[ ${ZISA_HAS_MPI} -eq 0 ]]
then
    zisa_dependencies=("ZisaCore" "ZisaMemory" "ZisaSFC" "ZisaTimeStepping")
else
    zisa_dependencies=("ZisaCore" "ZisaMemory" "ZisaSFC" "ZisaMPI" "ZisaTimeStepping")
fi

zisa_memory_root=$(realpath $(dirname $(readlink -f $0))/..)

install_dir=$(${zisa_memory_root}/bin/install_dir.sh $1 $2 --zisa_has_mpi=${ZISA_HAS_MPI})
source_dir=${install_dir}/sources
conan_file=${zisa_memory_root}/conanfile.txt

if [[ -f $conan_file ]]
then
   mkdir -p ${install_dir}/conan && cd ${install_dir}/conan
   conan install $conan_file -s compiler.libcxx=libstdc++11
fi

mkdir -p ${source_dir}
for dep in "${zisa_dependencies[@]}"
do
    src_dir=${source_dir}/$dep
    repo_url=git@github.com:1uc/${dep}.git

    if [[ -d ${src_dir} ]]
    then
        cd ${src_dir}

        if [[ -z $(git remote -v 2>/dev/null | grep ${repo_url}) ]]
        then
            echo "Failed to install ${dep} to ${src_dir}"
            exit -1

        else
            cd -
            rm -rf ${src_dir}
        fi
    fi

    git clone ${repo_url} ${src_dir}

    mkdir -p ${src_dir}/build-dep
    cd ${src_dir}/build-dep

    cmake -DCMAKE_INSTALL_PREFIX=${install_dir}/zisa \
          -DCMAKE_PREFIX_PATH=${install_dir}/zisa/lib/cmake/zisa \
          -DCMAKE_MODULE_PATH=${install_dir}/conan \
          -DCMAKE_PROJECT_${dep}_INCLUDE=${install_dir}/conan/conan_paths.cmake \
          -DZISA_HAS_MPI=${ZISA_HAS_MPI} \
          -DCMAKE_BUILD_TYPE=Release \
          ..

    cmake --build --parallel $(nproc) .
    cmake --install .
done

echo "The dependencies were installed at"
echo "    ${install_dir}"
