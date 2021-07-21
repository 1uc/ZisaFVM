// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/mpi/io/parallel_load_snapshot.hpp>

namespace zisa {

ParallelLoadSnapshot::ParallelLoadSnapshot(
    std::shared_ptr<FNG> fng,
    std::shared_ptr<HDF5UnstructuredFileDimensions> file_dimensions)
    : super(std::move(fng)), file_dims(std::move(file_dimensions)) {}

std::unique_ptr<HDF5Reader>
ParallelLoadSnapshot::pick_reader(const std::string &file_name) {
  return std::make_unique<HDF5UnstructuredReader>(file_name, file_dims);
}

}