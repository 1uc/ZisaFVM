#ifndef DUMP_SNAPSHOT_H_2AOQL
#define DUMP_SNAPSHOT_H_2AOQL

#include <zisa/grid/grid.hpp>
#include <zisa/io/file_name_generator.hpp>
#include <zisa/io/hdf5_writer.hpp>
#include <zisa/io/visualization.hpp>

namespace zisa {

/// Write the prognostic and diagnostic variables to the hard-disk.
template <class Model>
class DumpSnapshot : public Visualization {
public:
  DumpSnapshot(const Model &model,
               const std::shared_ptr<FileNameGenerator> &fng,
               const AllVariablesDimensions &dims);

protected:
  virtual void
  do_visualization(const AllVariables &all_variables,
                   const SimulationClock &simulation_clock) override;

protected:
  virtual std::unique_ptr<HDF5Writer> pick_writer(const std::string &file_name);

private:
  Model model;
  std::shared_ptr<FileNameGenerator> file_name_generator;
};

template <class Model>
class SerialDumpSnapshot : public Visualization {
public:
  SerialDumpSnapshot(const Model &model,
                     const std::shared_ptr<FileNameGenerator> &fng,
                     const std::shared_ptr<Grid> &grid);

protected:
  virtual void
  do_visualization(const AllVariables &all_variables,
                   const SimulationClock &simulation_clock) override;

private:
  Model model;
  std::shared_ptr<FileNameGenerator> file_name_generator;
  std::shared_ptr<Grid> grid;
};

} // namespace zisa

#endif
