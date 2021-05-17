#include <zisa/mpi/io/hdf5_unstructured_file_dimensions.hpp>

namespace zisa {

HDF5UnstructuredFileDimensions::HDF5UnstructuredFileDimensions(
    hsize_t n_cells_global,
    hsize_t n_cells_local,
    std::vector<hsize_t> ids,
    const MPI_Comm &mpi_comm)
    : n_cells_global(n_cells_global),
      n_cells_local(n_cells_local),
      ids(std::move(ids)),
      mpi_comm(mpi_comm) {}

}
