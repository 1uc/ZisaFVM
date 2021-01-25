#! /usr/bin/env bash

BSUB_2="bsub -n 2 -R span[hosts=1] -R rusage[mem=2000] -W 24:00"
BSUB_16="bsub -n 16 -R span[hosts=1] -R rusage[mem=2000] -W 24:00"
BSUB_16="bsub -n 32 -R span[hosts=1] -R rusage[mem=3800] -W 24:00"
BSUB="bsub -n 64 -R span[hosts=1] -R rusage[mem=7500] -W 24:00"

cd ${ZISA}/build-omp && make -j$(nproc) domain-decomposition && cd -

# for n in 2 4 8
# do
#     ${BSUB_2} ${ZISA}/build-omp/domain-decomposition -n ${n} -k 2\
#             --grid grids/scaling-0/grid.msh.h5 -o grids/scaling-0/partitioned/${n}
# done

# for n in 2 4 8 16 32 64
# do
#     ${BSUB_16} ${ZISA}/build-omp/domain-decomposition -n ${n} -k 16 \
#             --grid grids/scaling-1/grid.msh.h5 -o grids/scaling-1/partitioned/${n}
# done

# for n in 8 16 32 64 128 256 512
# do
#     ${BSUB} ${ZISA}/build-omp/domain-decomposition -n ${n} -k 128\
#             --grid grids/scaling-2/grid.msh.h5 -o grids/scaling-2/partitioned/${n}
# done

# for n in 32 64 128 256 512 1024 2048
# do
#     ${BSUB} ${ZISA}/build-omp/domain-decomposition -n ${n} -k 128 \
#             --grid grids/scaling-3/grid.msh.h5 -o grids/scaling-3/partitioned/${n}
# done


for n in 256
do
    ${BSUB_32} ${ZISA}/build-omp/domain-decomposition -n ${n} -k 32 \
            --grid grids/scaling-4/grid.msh.h5 -o grids/scaling-4/partitioned/${n}
done

# for n in 512 1024 2048
# do
#     ${BSUB} ${ZISA}/build-omp/domain-decomposition -n ${n} -k 64 \
#             --grid grids/scaling-4/grid.msh.h5 -o grids/scaling-4/partitioned/${n}
# done
