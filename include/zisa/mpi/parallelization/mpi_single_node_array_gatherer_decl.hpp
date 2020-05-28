#ifndef ZISA_MPI_SINGLE_NODE_ARRAY_GATHERER_DECL_HPP_KDWXN
#define ZISA_MPI_SINGLE_NODE_ARRAY_GATHERER_DECL_HPP_KDWXN

#include <zisa/config.hpp>

#include <zisa/mpi/mpi.hpp>
#include <zisa/mpi/parallelization/distributed_array_info.hpp>
#include <zisa/parallelization/array_gatherer.hpp>

namespace zisa {

template <class T, int n_dims>
class MPISingleNodeArrayGatherer : public ArrayGatherer<T, n_dims> {
private:
  using super = ArrayGatherer<T, n_dims>;

protected:
  using view_t = typename super::view_t;
  using const_view_t = typename super::const_view_t;

public:
  MPISingleNodeArrayGatherer(std::shared_ptr<DistributedArrayInfo> array_info,
                             MPI_Comm mpi_comm,
                             int mpi_tag);

  void send(const const_view_t &const_view) const override;
  void receive(const view_t &view) const override;
  void copy_local_patch(const view_t &global,
                        const const_view_t &local) const override;

  bool is_this_rank_gathering() const override;

private:
  std::shared_ptr<DistributedArrayInfo> array_info;

  int tag;
  int rank;
  int comm_size;
  int gather_rank = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
};

class MPISingleNodeArrayGathererFactory {
public:
  MPISingleNodeArrayGathererFactory(
      std::shared_ptr<DistributedArrayInfo> array_info,
      MPI_Comm comm,
      int mpi_base_tag)
      : array_info(std::move(array_info)),
        comm(comm),
        current_mpi_tag(mpi_base_tag) {}

  template <class T, int n_dims>
  MPISingleNodeArrayGatherer<T, n_dims> create_object() {
    return MPISingleNodeArrayGatherer<T, n_dims>(
        array_info, comm, current_mpi_tag++);
  }

  template <class T, int n_dims>
  std::unique_ptr<MPISingleNodeArrayGatherer<T, n_dims>> create_pointer() {
    return std::make_unique<MPISingleNodeArrayGatherer<T, n_dims>>(
        array_info, comm, current_mpi_tag++);
  }

private:
  std::shared_ptr<DistributedArrayInfo> array_info;
  MPI_Comm comm;
  int current_mpi_tag;
};

}

#endif // ZISA_MPI_SINGLE_NODE_ARRAY_GATHERER_DECL_HPP
