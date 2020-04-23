#include <zisa/mpi/parallelization/mpi_single_node_array_gatherer.hpp>

namespace zisa {

DistributedArrayInfo::DistributedArrayInfo(array<int_t, 1> partition)
    : partition(std::move(partition)) {}

}