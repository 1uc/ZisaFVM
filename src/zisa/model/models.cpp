#include <zisa/model/models.hpp>

namespace zisa {
void save_state(HDF5Writer &writer,
                const AllVariables &u,
                double t,
                int_t n_steps,
                const std::vector<std::string> &labels) {
  writer.write_scalar(t, "time");
  writer.write_scalar(n_steps, "n_steps");

  save(writer, u, labels);
}

std::pair<double, int_t> load_state(HDF5Reader &reader,
                                    AllVariables &u,
                                    const std::vector<std::string> &labels) {

  auto time = reader.read_scalar<double>("time");
  auto n_steps = reader.read_scalar<int_t>("n_steps");

  AllVariables::load(reader, u, labels);
  return {time, n_steps};
}

}
