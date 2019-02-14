/* Implementation for visualizing the flow.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-09-14
 */
#ifndef VISUALIZATION_INC_F9TKQD6Z
#define VISUALIZATION_INC_F9TKQD6Z

#include <zisa/io/tri_plot.hpp>
#include <zisa/io/visualization.hpp>
#include <zisa/model/all_variables.hpp>

namespace zisa {

template <class Model>
DumpSnapshot<Model>::DumpSnapshot(
    const Model &model,
    const std::shared_ptr<FileNameGenerator> &file_name_generator)
    : model(model), file_name_generator(file_name_generator) {}

template <class Model>
void DumpSnapshot<Model>::do_visualization(
    const AllVariables &all_variables,
    const SimulationClock &simulation_clock) {

  double t = simulation_clock.current_time();
  int n_steps = simulation_clock.current_step();

  auto writer = pick_writer(file_name_generator->next_name());

  save_state(*writer, model, *all_variables, t, n_steps);
}

template <class Model>
std::unique_ptr<HDF5Writer>
SerialDumpSnapshot<Model>::pick_writer(const std::string &file_name) {
  return std::make_unique<HDF5SerialWriter>(file_name);
}

} // namespace zisa
#endif /* end of include guard */
