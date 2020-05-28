#include <zisa/mpi/parallelization/distributed_array_info.hpp>

namespace zisa {

DistributedArrayInfo::DistributedArrayInfo(array<int_t, 1> partition)
    : partition(std::move(partition)) {}

std::shared_ptr<DistributedArrayInfo>
make_distributed_array_info(array<int_t, 1> partition) {
  return std::make_shared<DistributedArrayInfo>(std::move(partition));
}
}
