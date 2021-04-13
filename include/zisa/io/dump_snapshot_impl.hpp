/* Implementation for visualizing the flow.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-09-14
 */
#ifndef VISUALIZATION_INC_F9TKQD6Z
#define VISUALIZATION_INC_F9TKQD6Z
#include "dump_snapshot_decl.hpp"

#include <zisa/io/file_name_generator.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/grid_variables_impl.hpp>
#include <zisa/model/save_full_state.hpp>

namespace zisa {

template <class EOS>
DumpSnapshot<EOS>::DumpSnapshot(std::shared_ptr<LocalEOSState<EOS>> local_eos,
                                std::shared_ptr<FNG> fng)
    : local_eos(std::move(local_eos)), fng(std::move(fng)) {}

template <class EOS>
void DumpSnapshot<EOS>::do_visualization(
    const AllVariables &all_variables,
    const SimulationClock &simulation_clock) {

  auto t = simulation_clock.current_time();
  auto n_steps = simulation_clock.current_step();

  auto writer = pick_writer(fng->next_name());
  local_eos->compute(all_variables);
  save_full_state(*writer, *local_eos, all_variables, t, n_steps);
}

template <class EOS>
void DumpSnapshot<EOS>::do_steady_state(const AllVariables &steady_state) {
  auto writer = pick_writer(fng->steady_state());
  save(*writer, steady_state, all_labels<typename EOS::cvars_t>());
}

template <class EOS>
std::unique_ptr<HierarchicalWriter>
SerialDumpSnapshot<EOS>::pick_writer(const std::string &file_name) {
  return std::make_unique<HDF5SerialWriter>(file_name);
}

} // namespace zisa
#endif /* end of include guard */
