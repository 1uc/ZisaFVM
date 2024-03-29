// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_HDF5_UNSTRUCTURED_WRITER_HPP
#define ZISA_HDF5_UNSTRUCTURED_WRITER_HPP

#include <zisa/config.hpp>

#include <numeric>
#include <zisa/io/hdf5_writer.hpp>
#include <zisa/mpi/io/hdf5_unstructured_file_dimensions.hpp>
#include <zisa/mpi/mpi.hpp>
#include <zisa/utils/integer_cast.hpp>

namespace zisa {

/// `n_cells` is the number of cells in the patch, including ghost-cells.
std::shared_ptr<HDF5UnstructuredFileDimensions>
make_hdf5_unstructured_file_dimensions(int_t n_cells,
                                       std::vector<hsize_t> ids,
                                       const MPI_Comm &mpi_comm);

enum class HDF5Access { overwrite, modify };

class HDF5ParallelWriter : public HDF5Writer {
public:
protected:
  void do_write_scalar(void const *addr,
                       const HDF5DataType &data_type,
                       const std::string &tag) override;

  void do_write_string(const std::string &data,
                       const std::string &tag) override;

  virtual bool is_serial_writer() const = 0;
};

class HDF5UnstructuredWriter : public HDF5ParallelWriter {
public:
  HDF5UnstructuredWriter(
      const std::string &filename,
      std::shared_ptr<HDF5UnstructuredFileDimensions> file_dims,
      const HDF5Access &access = HDF5Access::overwrite);

protected:
  void do_write_array(void const *data,
                      const HDF5DataType &data_type,
                      const std::string &tag,
                      int rank,
                      hsize_t const *dims) override;

  bool is_serial_writer() const override;

private:
  std::vector<hsize_t> local_count(int rank_, hsize_t const *dims) const;
  std::vector<hsize_t> global_ids(int rank_, hsize_t const *dims) const;

private:
  std::shared_ptr<HDF5UnstructuredFileDimensions> file_dims;
};

/// Read data from HDF5 file sequentially.
class HDF5UnstructuredReader : public HDF5Reader {
private:
  using super = HDF5Reader;

public:
  explicit HDF5UnstructuredReader(
      const std::string &filename,
      std::shared_ptr<HDF5UnstructuredFileDimensions> file_dims_);

  virtual std::vector<hsize_t>
  do_hdf5_dims(const std::string &tag) const override;

  virtual void do_read_array(void *data,
                             const HDF5DataType &data_type,
                             const std::string &tag) const override;

  virtual void do_read_scalar(void *data,
                              const HDF5DataType &data_type,
                              const std::string &tag) const override;

  virtual std::string do_read_string(const std::string &tag) const override;

private:
  std::shared_ptr<HDF5UnstructuredFileDimensions> file_dims;
};

}
#endif // ZISA_HDF5_UNSTRUCTURED_WRITER_HPP
