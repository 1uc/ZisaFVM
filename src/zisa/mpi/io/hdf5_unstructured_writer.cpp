#include <zisa/mpi/io/hdf5_unstructured_writer.hpp>

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

HDF5UnstructuredWriter::HDF5UnstructuredWriter(
    const std::string &filename,
    std::shared_ptr<HDF5UnstructuredFileDimensions> file_dims_)
    : file_dims(std::move(file_dims_)) {
  auto lock = std::lock_guard(hdf5_mutex);

  auto h5_plist = zisa::H5P::create(H5P_FILE_ACCESS);
  zisa::H5P::set_fapl_mpio(h5_plist, file_dims->mpi_comm, MPI_INFO_NULL);

  hid_t h5_file = zisa::H5F::create(
      filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, h5_plist);
  zisa::H5P::close(h5_plist);

  file.push(h5_file);
}

void HDF5UnstructuredWriter::do_write_scalar(const void *addr,
                                             const HDF5DataType &data_type,
                                             const std::string &tag) const {
  // create a scalar data space.
  hid_t dataspace = zisa::H5S::create(H5S_SCALAR);

  // Create dataspace, everyone needs to participate here! This way the
  // meta-data stays consistent among all processes.
  hid_t dataset = zisa::H5D::create(file.top(),
                                    tag.c_str(),
                                    data_type(),
                                    dataspace,
                                    H5P_DEFAULT,
                                    H5P_DEFAULT,
                                    H5P_DEFAULT);

  if (zisa::mpi::rank(file_dims->mpi_comm) == 0) {
    // we want only one process to write, hence we need `independent` mode.
    hid_t h5_properties = zisa::H5P::create(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(h5_properties, H5FD_MPIO_INDEPENDENT);

    // write the scalar, but only by one process
    zisa::H5D::write(
        dataset, data_type(), H5S_ALL, H5S_ALL, h5_properties, addr);

    zisa::H5P::close(h5_properties);
  }

  // close
  zisa::H5D::close(dataset);
  zisa::H5S::close(dataspace);

  MPI_Barrier(file_dims->mpi_comm);
}

void HDF5UnstructuredWriter::do_write_string(const std::string &data,
                                             const std::string &tag) const {
  // strings can be stored as 1d-arrays of characters.
  // don't forget the null-character at the end of 'data.c_str()'.
  hsize_t dims[1] = {data.size() + 1};
  hid_t dataspace = H5Screate_simple(1, dims, nullptr);

  // this type of characters
  HDF5DataType data_type = make_hdf5_data_type<char>();

  // Create dataspace, everyone needs to participate here! This way the
  // meta-data stays consistent among all processes.
  hid_t dataset = zisa::H5D::create(file.top(),
                                    tag.c_str(),
                                    data_type(),
                                    dataspace,
                                    H5P_DEFAULT,
                                    H5P_DEFAULT,
                                    H5P_DEFAULT);

  if (zisa::mpi::rank(file_dims->mpi_comm) == 0) {
    // we want only one process to write, hence we need `independent` mode.
    hid_t h5_properties = zisa::H5P::create(H5P_DATASET_XFER);
    zisa::H5P::set_dxpl_mpio(h5_properties, H5FD_MPIO_INDEPENDENT);

    // write the string, but only by one process
    zisa::H5D::write(
        dataset, data_type(), H5S_ALL, H5S_ALL, h5_properties, data.c_str());

    zisa::H5P::close(h5_properties);
  }

  // close
  zisa::H5D::close(dataset);
  zisa::H5S::close(dataspace);
}

void HDF5UnstructuredWriter::do_write_array(const void *data,
                                            const HDF5DataType &data_type,
                                            const std::string &tag,
                                            int rank,
                                            hsize_t const *dims) const {

  // assert incoming shape.

  auto global_dims = std::vector<hsize_t>(integer_cast<size_t>(rank));
  global_dims[0] = file_dims->n_cells_global;
  for (size_t r = 1; r < integer_cast<size_t>(rank); ++r) {
    global_dims[r] = dims[r];
  }

  // create filespace
  hid_t h5_filespace
      = zisa::H5S::create_simple(rank, global_dims.data(), nullptr);

  // create a property list
  hid_t h5_plist = zisa::H5P::create(H5P_DATASET_CREATE);

  // create a dataset
  hid_t h5_dataset = zisa::H5D::create(file.top(),
                                       tag.c_str(),
                                       H5T_NATIVE_DOUBLE,
                                       h5_filespace,
                                       H5P_DEFAULT,
                                       h5_plist,
                                       H5P_DEFAULT);

  zisa::H5P::close(h5_plist);
  zisa::H5S::close(h5_filespace);

  // NOTE: HDF5 tutorial suggest this extra closing and reopening step of the
  // file space
  h5_filespace = zisa::H5D::get_space(h5_dataset);

  // select the corresponding slab on disk
  auto gids = global_ids(rank, dims);
  H5S::select_elements(h5_filespace,
                       H5S_SELECT_SET,
                       gids.size() / integer_cast<size_t>(rank),
                       gids.data());

  // select only part of the memspace
  hid_t h5_memspace = zisa::H5S::create_simple(rank, dims, nullptr);

  auto offset = std::vector<hsize_t>(integer_cast<size_t>(rank), hsize_t(0));
  auto count = local_count(rank, dims);
  zisa::H5S::select_hyperslab(h5_memspace,
                              H5S_SELECT_SET,
                              offset.data(),
                              nullptr,
                              count.data(),
                              nullptr);

  // property list for collective MPI write
  h5_plist = zisa::H5P::create(H5P_DATASET_XFER);
  zisa::H5P::set_dxpl_mpio(h5_plist, H5FD_MPIO_COLLECTIVE);

  // write slab to file
  zisa::H5D::write(
      h5_dataset, data_type(), h5_memspace, h5_filespace, h5_plist, data);

  // close everything
  zisa::H5S::close(h5_filespace);
  zisa::H5S::close(h5_memspace);
  zisa::H5P::close(h5_plist);
  zisa::H5D::close(h5_dataset);
}

std::vector<hsize_t>
HDF5UnstructuredWriter::local_count(int rank_,
                                    hsize_t const *const dims) const {
  auto rank = integer_cast<size_t>(rank_);

  std::vector<hsize_t> count(rank);
  std::copy(dims, dims + rank, count.data());
  count[0] = file_dims->ids.size();

  return count;
}

hsize_t trailing_product(size_t rank, hsize_t const *const dims) {
  return std::accumulate(dims + 1,
                         dims + rank,
                         hsize_t(1),
                         [](hsize_t a, hsize_t b) { return a * b; });
}

std::vector<hsize_t> global_ids(const std::vector<hsize_t> &ids,
                                size_t rank,
                                hsize_t const *const dims) {

  auto nids = ids.size();
  auto nflat = trailing_product(rank, dims);
  auto gids = std::vector<hsize_t>(nids * nflat * rank);

  auto ld = std::vector<hsize_t>(rank);
  ld.back() = 1;
  for (size_t r = rank - 1; r > 0; --r) {
    ld[r - 1] = dims[r] * ld[r];
  }

  for (hsize_t i = 0; i < nids; ++i) {
    for (hsize_t j = 0; j < nflat; ++j) {
      gids[(i * ld[0] + j) * rank + 0] = ids[i];
      for (size_t r = 1; r < rank; ++r) {
        gids[(i * ld[0] + j) * rank + r] = (j % ld[r - 1]) / ld[r];
      }
    }
  }

  return gids;
}

std::vector<hsize_t>
HDF5UnstructuredWriter::global_ids(int rank, hsize_t const *const dims) const {
  return zisa::global_ids(file_dims->ids, integer_cast<size_t>(rank), dims);
}

std::shared_ptr<HDF5UnstructuredFileDimensions>
make_hdf5_unstructured_file_dimensions(int_t n_cells,
                                       std::vector<hsize_t> ids,
                                       const MPI_Comm &mpi_comm) {
  size_t n_local_cells = ids.size();
  size_t n_total_cells = 0;
  MPI_Allreduce(
      &n_local_cells, &n_total_cells, 1, MPI_SIZE_T, MPI_SUM, mpi_comm);

  return std::make_shared<HDF5UnstructuredFileDimensions>(
      integer_cast<hsize_t>(n_total_cells),
      integer_cast<hsize_t>(n_cells),
      std::move(ids),
      mpi_comm);
}

HDF5UnstructuredReader::HDF5UnstructuredReader(
    const std::string &filename,
    std::shared_ptr<HDF5UnstructuredFileDimensions> file_dims_)
    : file_dims(std::move(file_dims_)) {

  // define a properties list for a parallel HDF5 file.
  hid_t h5_plist = zisa::H5P::create(H5P_FILE_ACCESS);
  zisa::H5P::set_fapl_mpio(h5_plist, file_dims->mpi_comm, MPI_INFO_NULL);

  hid_t h5_file = zisa::H5F::open(filename.c_str(), H5F_ACC_RDONLY, h5_plist);
  zisa::H5P::close(h5_plist);

  file.push(h5_file);
  path.push_back(filename);
}

std::vector<hsize_t>
HDF5UnstructuredReader::do_dims(const std::string &tag) const {
  // TODO this is duplicated in `HDF5SerialReader`.
  hid_t dataset = open_dataset(tag);
  hid_t dataspace = get_dataspace(dataset);

  auto rank = static_cast<int_t>(H5S::get_simple_extent_ndims(dataspace));
  std::vector<hsize_t> dims(rank);

  H5S::get_simple_extent_dims(dataspace, &(dims[0]), nullptr);
  dims[0] = file_dims->n_cells_local;

  zisa::H5D::close(dataset);
  zisa::H5S::close(dataspace);

  return dims;
}

void HDF5UnstructuredReader::do_read_array(void *data,
                                           const HDF5DataType &data_type,
                                           const std::string &tag) const {

  auto local_dims = dims(tag);
  auto rank = local_dims.size();
  local_dims[0] = file_dims->ids.size();

  hid_t dataset = open_dataset(tag);
  hid_t filespace = get_dataspace(dataset);

  hid_t memspace = zisa::H5S::create_simple(
      integer_cast<int>(rank), local_dims.data(), nullptr);

  hid_t plist = zisa::H5P::create(H5P_DATASET_XFER);
  zisa::H5P::set_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

  auto gids = global_ids(file_dims->ids, rank, local_dims.data());
  H5S::select_elements(
      filespace, H5S_SELECT_SET, gids.size() / rank, gids.data());

  // read slab from file
  zisa::H5D::read(dataset, data_type(), memspace, filespace, H5P_DEFAULT, data);

  zisa::H5P::close(plist);
  zisa::H5S::close(filespace);
  zisa::H5S::close(memspace);
  zisa::H5D::close(dataset);
}

void HDF5UnstructuredReader::do_read_scalar(void *data,
                                            const HDF5DataType &data_type,
                                            const std::string &tag) const {

  auto dataset = open_dataset(tag);
  auto dataspace = get_dataspace(dataset);

  auto properties = zisa::H5P::create(H5P_DATASET_XFER);
  zisa::H5P::set_dxpl_mpio(properties, H5FD_MPIO_COLLECTIVE);

  zisa::H5D::read(dataset, data_type(), H5S_ALL, H5S_ALL, properties, data);

  zisa::H5D::close(dataset);
  zisa::H5S::close(dataspace);
  zisa::H5P::close(properties);
}

std::string HDF5UnstructuredReader::do_read_string(const std::string &) const {
  LOG_ERR("Implement first.");
}
}
