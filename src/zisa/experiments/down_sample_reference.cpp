#include <zisa/experiments/down_sample_reference.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/io/file_manipulation.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/math/reference_solution.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/models.hpp>

namespace zisa {

void down_sample_euler_reference(
    const ReferenceSolution &reference_solution,
    const std::vector<std::string> &coarse_grid_paths,
    const std::string &filename) {

  for (const auto &grid_name : coarse_grid_paths) {
    auto coarse_grid = load_grid(grid_name, /* quad_deg = */ 4);
    auto all_vars_coarse = reference_solution.average(*coarse_grid);

    std::string stem = zisa::stem(zisa::basename(grid_name));
    std::string output_name
        = string_format("down_sampled/%s/%s", stem.c_str(), filename.c_str());

    auto writer = HDF5SerialWriter(output_name);
    save(writer, *all_vars_coarse, all_labels<euler_var_t>());
  }
}

}