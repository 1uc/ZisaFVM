#ifndef ZISA_HDF5_UNSTRUCTURED_WRITER_HPP
#define ZISA_HDF5_UNSTRUCTURED_WRITER_HPP

#include <zisa/config.hpp>

#include <numeric>
#include <zisa/io/hdf5_writer.hpp>
#include <zisa/mpi/mpi.hpp>
#include <zisa/utils/integer_cast.hpp>

namespace zisa {

struct HDF5UnstructuredFileDimensions {
  hsize_t n_cells;
  std::vector<hsize_t> ids;
  MPI_Comm mpi_comm;

  HDF5UnstructuredFileDimensions(hsize_t n_cells,
                                 std::vector<hsize_t> ids,
                                 const MPI_Comm &mpi_comm);
};

std::shared_ptr<HDF5UnstructuredFileDimensions>
make_hdf5_unstructured_file_dimensions(std::vector<hsize_t> ids,
                                       const MPI_Comm &mpi_comm);

class HDF5UnstructuredWriter : public HDF5Writer {
public:
  HDF5UnstructuredWriter(
      const std::string &filename,
      std::shared_ptr<HDF5UnstructuredFileDimensions> file_dims);

protected:
  void do_write_array(void const *data,
                      const HDF5DataType &data_type,
                      const std::string &tag,
                      int rank,
                      hsize_t const *dims) const override {

    // assert incoming shape.

    auto global_dims = std::vector<hsize_t>(integer_cast<size_t>(rank));
    global_dims[0] = n_cells;
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

  void do_write_scalar(void const *addr,
                       const HDF5DataType &data_type,
                       const std::string &tag) const override {
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

    if (zisa::mpi::rank(mpi_comm) == 0) {
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

    MPI_Barrier(mpi_comm);
  }

  void do_write_string(const std::string &data,
                       const std::string &tag) const override {
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

    if (zisa::mpi::rank(mpi_comm) == 0) {
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

private:
  std::vector<hsize_t> local_count(int rank_, hsize_t const *const dims) const {
    auto rank = integer_cast<size_t>(rank_);

    std::vector<hsize_t> count(rank);
    std::copy(dims, dims + rank, count.data());
    count[0] = ids.size();

    return count;
  }

  std::vector<hsize_t> global_ids(int rank_, hsize_t const *const dims) const {
    auto rank = integer_cast<size_t>(rank_);

    auto nids = ids.size();
    auto nflat = trailing_product(rank_, dims);
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

  hsize_t trailing_product(int rank, hsize_t const *const dims) const {
    return std::accumulate(dims + 1,
                           dims + rank,
                           hsize_t(1),
                           [](hsize_t a, hsize_t b) { return a * b; });
  }

private:
  std::shared_ptr<HDF5UnstructuredFileDimensions> file_dims;

  hsize_t n_cells;
  const std::vector<hsize_t> &ids;
  MPI_Comm mpi_comm;
};
}
#endif // ZISA_HDF5_UNSTRUCTURED_WRITER_HPP
