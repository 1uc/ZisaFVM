#include <zisa/mpi/io/hdf5_unstructured_writer.hpp>

namespace zisa {

HDF5UnstructuredFileDimensions::HDF5UnstructuredFileDimensions(
    hsize_t n_cells, std::vector<hsize_t> ids, const MPI_Comm &mpi_comm)
    : n_cells(n_cells), ids(std::move(ids)), mpi_comm(mpi_comm) {}

HDF5UnstructuredWriter::HDF5UnstructuredWriter(
    const std::string &filename,
    std::shared_ptr<HDF5UnstructuredFileDimensions> file_dims)
    : file_dims(std::move(file_dims)),
      n_cells(this->file_dims->n_cells),
      ids(this->file_dims->ids),
      mpi_comm(this->file_dims->mpi_comm) {
  auto lock = std::lock_guard(hdf5_mutex);

  auto h5_plist = zisa::H5P::create(H5P_FILE_ACCESS);
  zisa::H5P::set_fapl_mpio(h5_plist, mpi_comm, MPI_INFO_NULL);

  hid_t h5_file = zisa::H5F::create(
      filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, h5_plist);
  zisa::H5P::close(h5_plist);

  file.push(h5_file);
}

std::shared_ptr<HDF5UnstructuredFileDimensions>
make_hdf5_unstructured_file_dimensions(std::vector<hsize_t> ids,
                                       MPI_Comm const &mpi_comm) {
  size_t n_local_cells = ids.size();
  size_t n_total_cells = 0;
  MPI_Allreduce(
      &n_local_cells, &n_total_cells, 1, MPI_SIZE_T, MPI_SUM, mpi_comm);

  return std::make_shared<HDF5UnstructuredFileDimensions>(
      integer_cast<hsize_t>(n_total_cells), std::move(ids), mpi_comm);
}

}
