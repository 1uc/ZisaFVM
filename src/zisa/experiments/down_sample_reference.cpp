#include <zisa/experiments/down_sample_reference.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/models.hpp>

namespace zisa {
void down_sample_reference(const InputParameters &input_params) {

  const auto &params = input_params["down-sample"];

  auto fine_grid = load_gmsh(params["reference"]["grid"]);

  auto reader = HDF5SerialReader(params["reference"]["data"]);

  std::string model = reader.read_string("model");
  LOG_ERR_IF(model != "Euler",
             string_format("Unknown Model. [%s]", model.c_str()));

  auto labels = all_labels<euler_var_t>();
  auto fine_data = AllVariables::load(reader, labels);
}
}