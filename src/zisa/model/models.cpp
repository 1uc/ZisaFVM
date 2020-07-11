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
}
