/* Functionality shared by all Models.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2015-01-20
 */
#ifndef MODELS_H_KL3BI6CZ
#define MODELS_H_KL3BI6CZ

#include <string>
#include <vector>

#include <zisa/config.hpp>
#include <zisa/model/all_variables.hpp>

namespace zisa {
template <class cvars_t> std::vector<std::string> all_labels(void) {
  std::vector<std::string> ret(cvars_t::size());
  for (int k = 0; k < cvars_t::size(); ++k) {
    ret[k] = cvars_t::labels(k);
  }

  return ret;
}

// template <class Model>
// void save_state(HDF5Writer &writer, const Model &model, const AllVariables &u,
//                 double t, int n_steps) {
//   model.save_parameters(writer);

//   writer.write_scalar(t, "time");
//   writer.write_scalar(n_steps, "n_steps");

//   u.save(writer, all_labels<typename Model::cvars_t>());
// }
} // namespace tyr
#endif /* end of include guard: MODELS_H_KL3BI6CZ */
