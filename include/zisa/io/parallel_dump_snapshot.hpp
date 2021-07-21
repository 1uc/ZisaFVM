// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_PARALLEL_DUMP_SNAPSHOT_HPP_UUYBW
#define ZISA_PARALLEL_DUMP_SNAPSHOT_HPP_UUYBW

#include <zisa/io/dump_snapshot.hpp>
#include <zisa/mpi/io/hdf5_unstructured_writer.hpp>

namespace zisa {

template <class EOS>
class ParallelDumpSnapshot : public DumpSnapshot<EOS> {
private:
  using super = DumpSnapshot<EOS>;

public:
  ParallelDumpSnapshot(
      std::shared_ptr<LocalEOSState<EOS>> local_eos,
      std::shared_ptr<FNG> fng,
      std::shared_ptr<HDF5UnstructuredFileDimensions> file_dimensions)
      : super(std::move(local_eos), std::move(fng)),
        file_dims(std::move(file_dimensions)) {}

protected:
  virtual std::unique_ptr<HierarchicalWriter>
  pick_writer(const std::string &file_name) override {
    return std::make_unique<HDF5UnstructuredWriter>(file_name, file_dims);
  }

private:
  std::shared_ptr<HDF5UnstructuredFileDimensions> file_dims;
};

}
#endif
