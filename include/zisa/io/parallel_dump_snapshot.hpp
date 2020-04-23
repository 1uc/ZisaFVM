#ifndef ZISA_PARALLEL_DUMP_SNAPSHOT_HPP_UUYBW
#define ZISA_PARALLEL_DUMP_SNAPSHOT_HPP_UUYBW

#include <zisa/io/dump_snapshot.hpp>
#include <zisa/mpi/io/hdf5_unstructured_writer.hpp>

namespace zisa {

template <class Model>
class ParallelDumpSnapshot : public DumpSnapshot<Model> {
private:
  using super = DumpSnapshot<Model>;

public:
  ParallelDumpSnapshot(
      std::shared_ptr<Model> model,
      std::shared_ptr<FileNameGenerator> file_name_generator,
      std::shared_ptr<HDF5UnstructuredFileDimensions> file_dimensions)
      : super(std::move(model), std::move(file_name_generator)),
        file_dims(std::move(file_dimensions)) {}

protected:
  virtual std::unique_ptr<HDF5Writer>
  pick_writer(const std::string &file_name) override {
    return std::make_unique<HDF5UnstructuredWriter>(file_name, file_dims);
  }

private:
  std::shared_ptr<HDF5UnstructuredFileDimensions> file_dims;
};

}
#endif
