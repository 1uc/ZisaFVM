// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef DUMP_SNAPSHOT_H_2AOQL
#define DUMP_SNAPSHOT_H_2AOQL

#include <zisa/grid/grid.hpp>
#include <zisa/io/file_name_generator.hpp>
#include <zisa/io/hdf5_writer.hpp>
#include <zisa/io/visualization.hpp>
#include <zisa/model/all_variables_fwd.hpp>
#include <zisa/model/grid_variables_decl.hpp>
#include <zisa/model/local_eos_state.hpp>

namespace zisa {

/// Write the prognostic and diagnostic variables to the hard-disk.
template <class EOS>
class DumpSnapshot : public Visualization {
public:
  DumpSnapshot(std::shared_ptr<LocalEOSState<EOS>> eos,
               std::shared_ptr<FNG> fng);

protected:
  virtual void
  do_visualization(const AllVariables &all_variables,
                   const SimulationClock &simulation_clock) override;

  virtual void do_steady_state(const AllVariables &steady_state) override;

protected:
  virtual std::unique_ptr<HierarchicalWriter> pick_writer(const std::string &file_name)
      = 0;

private:
  std::shared_ptr<LocalEOSState<EOS>> local_eos;
  std::shared_ptr<FNG> fng;
};

template <class EOS>
class SerialDumpSnapshot : public DumpSnapshot<EOS> {
private:
  using super = DumpSnapshot<EOS>;

public:
  using super::super;

protected:
  virtual std::unique_ptr<HierarchicalWriter>
  pick_writer(const std::string &file_name) override;
};

} // namespace zisa

#endif
