// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_PARALLEL_LOAD_SNAPSHOT_HPP_BUEWI
#define ZISA_PARALLEL_LOAD_SNAPSHOT_HPP_BUEWI

#include <zisa/config.hpp>

#include <zisa/io/file_name_generator.hpp>
#include <zisa/io/load_snapshot.hpp>
#include <zisa/mpi/io/hdf5_unstructured_writer.hpp>

namespace zisa {

class ParallelLoadSnapshot : public LoadSnapshot {
private:
  using super = LoadSnapshot;

public:
  ParallelLoadSnapshot(
      std::shared_ptr<FNG> fng,
      std::shared_ptr<HDF5UnstructuredFileDimensions> file_dimensions);

protected:
  virtual std::unique_ptr<HDF5Reader>
  pick_reader(const std::string &file_name) override;

private:
  std::shared_ptr<HDF5UnstructuredFileDimensions> file_dims;
};

}

#endif // ZISA_PARALLEL_LOAD_SNAPSHOT_HPP
