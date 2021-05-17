#ifndef HDF5_UNSTRUCTURED_FILE_DIMENSIONS_HPP_3KGKNX27
#define HDF5_UNSTRUCTURED_FILE_DIMENSIONS_HPP_3KGKNX27

#include <zisa/io/hdf5.hpp>
#include <zisa/mpi/mpi.hpp>

namespace zisa {

struct HDF5UnstructuredFileDimensions {
  hsize_t n_cells_global;
  hsize_t n_cells_local;
  std::vector<hsize_t> ids;
  MPI_Comm mpi_comm;

  HDF5UnstructuredFileDimensions(hsize_t n_cells_global,
                                 hsize_t n_cells_local,
                                 std::vector<hsize_t> ids,
                                 const MPI_Comm &mpi_comm);
};

}
#endif
