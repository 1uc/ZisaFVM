#ifndef ZISA_MPI_SINGLE_NODE_ARRAY_SCATTERER_DECL_HPP_XCJQPO
#define ZISA_MPI_SINGLE_NODE_ARRAY_SCATTERER_DECL_HPP_XCJQPO

#include <zisa/config.hpp>

#include <zisa/mpi/parallelization/distributed_array_info.hpp>
#include <zisa/parallelization/array_scatterer.hpp>

namespace zisa {

template <class T, int n_dims>
class MPISingleNodeArrayScatterer : public ArrayScatterer<T, n_dims> {
private:
  using super = ArrayScatterer<T, n_dims>;

protected:
  using view_t = typename super::view_t;
  using const_view_t = typename super::const_view_t;

public:
  MPISingleNodeArrayScatterer(std::shared_ptr<DistributedArrayInfo> array_info,
                              MPI_Comm mpi_comm,
                              int mpi_tag);

  void send(const const_view_t &const_view) const override;
  void receive(const view_t &view) const override;
  void copy_local_patch(const view_t &local,
                        const const_view_t &global) const override;

  bool is_this_rank_scattering() const override;

private:
  std::shared_ptr<DistributedArrayInfo> array_info;

  int tag;
  int rank;
  int comm_size;
  int scatter_rank = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
};
}

#endif // ZISA_MPI_SINGLE_NODE_ARRAY_SCATTERER_DECL_HPP
