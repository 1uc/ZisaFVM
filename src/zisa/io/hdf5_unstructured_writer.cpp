#include <zisa/io/hdf5_unstructured_writer.hpp>

namespace zisa {

HDF5UnstructuredWriter::HDF5UnstructuredWriter(
    const std::string &filename,
    HDF5UnstructuredFileDimensions file_dims,
    const MPI_Comm &mpi_comm)
    : n_cells(file_dims.n_cells),
      ids(std::move(file_dims.ids)),
      mpi_comm(mpi_comm) {

  auto lock = std::lock_guard(hdf5_mutex);

  auto h5_plist = zisa::H5P::create(H5P_FILE_ACCESS);
  zisa::H5P::set_fapl_mpio(h5_plist, mpi_comm, MPI_INFO_NULL);

  hid_t h5_file = zisa::H5F::create(
      filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, h5_plist);
  zisa::H5P::close(h5_plist);

  file.push(h5_file);
}

}