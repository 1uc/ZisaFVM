#ifndef ZISA_GATHERED_VIS_INFO_HPP_CDUIEW
#define ZISA_GATHERED_VIS_INFO_HPP_CDUIEW

#include <zisa/config.hpp>

#include <memory>

#include <zisa/math/permutation.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/mpi/mpi.hpp>

namespace zisa {

/// Datastructure which describes Multiwriter data layout.
/** The idea is to collect data the on a fixed number of processing elements /
 *  processes / MPI-tasks (PEs). Then use HDF5UnstructuredWriter to write the
 *  data to disk.
 *
 *  A 'visualizing PE' is one that will participate in the writing of the HDF5
 *  file.
 *
 *  `vis_boundaries` this is the cummulative sum of the number of cell each
 *      visualizing PE owns.
 *
 *  `vis_file_ids` these are the global IDs of the cell. This tells
 *   HDF5UnstructuredWriter where to store the gathered chunk of data. They are
 *   sorted in ascending order.
 *
 *  `permutation` this permutation must be applied to the gathered data to
 *   ensure it's in the same order as `vis_file_ids`.
 *
 *   `n_local_cells` number of cell on this PE that need to be sent away to be
 *    written to disk.
 *
 *   `vis_comm` on these MPI comms one will gather the data.
 *   `h5_comm` MPI comm spanning all 'visualizing' PEs.
 * */
struct GatheredVisInfo {
  array<int_t, 1> vis_boundaries;
  array<int_t, 1> vis_file_ids;
  std::shared_ptr<Permutation> permutation;

  inline int_t n_vis_cells() const {
    if (h5_comm != MPI_COMM_NULL) {
      return vis_boundaries[vis_boundaries.size() - 1];
    } else {
      return 0;
    }
  }

  int_t n_local_cells;

  MPI_Comm vis_comm;
  MPI_Comm h5_comm;
};

}

#endif // ZISA_GATHERED_VIS_INFO_HPP
