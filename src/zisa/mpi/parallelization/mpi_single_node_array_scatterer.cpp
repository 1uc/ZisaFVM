#include <zisa/mpi/parallelization/mpi_single_node_array_scatterer.hpp>

namespace zisa {
MPISingleNodeArrayScattererFactory::MPISingleNodeArrayScattererFactory(
    std::shared_ptr<DistributedArrayInfo> array_info,
    MPI_Comm comm,
    int mpi_base_tag)
    : array_info(std::move(array_info)),
      comm(comm),
      current_mpi_tag(mpi_base_tag) {}

std::shared_ptr<MPISingleNodeArrayScattererFactory>
make_mpi_single_node_array_scatterer_factory(const GatheredVisInfo &vis_info) {
  auto vis_comm = vis_info.vis_comm;
  auto vis_tag = ZISA_MPI_TAG_SCATTERED_VIS;

  auto h5_comm = vis_info.h5_comm;

  if (h5_comm != MPI_COMM_NULL) {
    auto darray_info = make_distributed_array_info(vis_info.vis_boundaries);
    return std::make_shared<MPISingleNodeArrayScattererFactory>(
        darray_info, vis_comm, vis_tag + 1);
  } else {
    return std::make_shared<MPISingleNodeArrayScattererFactory>(
        nullptr, vis_comm, vis_tag + 1);
  }
}

}