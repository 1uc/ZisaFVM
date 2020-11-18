#include <zisa/model/load_full_state.hpp>

namespace zisa {

std::pair<double, int_t> load_full_state(HDF5Reader &reader,
                                         AllVariables &all_variables) {

  auto labels = all_labels<typename Euler::cvars_t>();
  return load_state(reader, all_variables, labels);
}

}