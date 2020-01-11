#ifndef ZISA_MPI_SINGLE_NODE_ARRAY_GATHERER_DECL_HPP_KDWXN
#define ZISA_MPI_SINGLE_NODE_ARRAY_GATHERER_DECL_HPP_KDWXN

#include <zisa/parallelization/array_gatherer.hpp>
#include <zisa/parallelization/mpi.hpp>

namespace zisa {

struct DistributedArrayInfo {
  explicit DistributedArrayInfo(array<int_t, 1> partition);

  // Rank `p` has all data with indices [partition[p], partition[p+1]).
  array<int_t, 1> partition;
};

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


}

#endif // ZISA_MPI_SINGLE_NODE_ARRAY_GATHERER_DECL_HPP
