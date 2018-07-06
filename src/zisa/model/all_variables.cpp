/* Data-structure for all cell-centered variables exposed to Tyr.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-06-13
 */
#ifndef ALL_VARIABLES_CPP_WLARDDKM
#define ALL_VARIABLES_CPP_WLARDDKM

#include <zisa/config.hpp>
#include <zisa/model/all_variables.hpp>

namespace zisa {
void GridVariables::split_save(HDF5Writer &writer,
                               const std::vector<std::string> &labels) const {
  array<double, 1> component({shape(0)}, device_type::cpu);
  for (int_t k = 0; k < labels.size(); ++k) {
    for (int_t i = 0; i < shape(0); ++i) {
      component[i] = (*this)(i, k);
    }

    zisa::save(writer, component, labels.at(k));
  }
}

// void GridVariables::split_load(HDF5Reader &reader,
//                                const std::vector<std::string> &labels) {
//   Array<double, 1> component;

//   for (int_t k = 0; k < labels.size(); ++k) {
//     component.load(reader, labels.at(k));

//     for (int_t i = 0; i < this->shape(1); ++i) {
//       (*this)(i, k) = component[i];
//     }
//   }
// }

AllVariables::AllVariables(const AllVariablesDimensions &dims,
                           const device_type &device) {
  allocate(dims, device);
}

void AllVariables::allocate(const AllVariablesDimensions &dims,
                            const device_type &device) {
  conserved_variables
      = GridVariables({dims.n_cells, dims.n_conserved_variables}, device);
  advected_variables
      = GridVariables({dims.n_cells, dims.n_advected_variables}, device);
}

double AllVariables::operator[](int_t i) const {
  int n_conserved_elements = conserved_variables.size();

  if (i < n_conserved_elements) {
    return conserved_variables[i];
  } else {
    return advected_variables[i - n_conserved_elements];
  }
}

double &AllVariables::operator[](int_t i) {
  int n_conserved_elements = conserved_variables.size();

  if (i < n_conserved_elements) {
    return conserved_variables[i];
  } else {
    return advected_variables[i - n_conserved_elements];
  }
}

int_t AllVariables::size() const {
  return conserved_variables.size() + advected_variables.size();
}

AllVariablesDimensions AllVariables::dims(void) const {
  AllVariablesDimensions dims;

  dims.n_cells = conserved_variables.shape(0);
  dims.n_conserved_variables = conserved_variables.shape(1);
  dims.n_advected_variables = advected_variables.shape(1);

  return dims;
}

void AllVariables::save(HDF5Writer &writer,
                        const std::vector<std::string> &labels) const {
  // conserved variables
  conserved_variables.split_save(writer, labels);

  int n_advected_variables = advected_variables.shape(1);
  assert(n_advected_variables == 0);
  // // advected variables
  // int n_advected_variables = advected_variables.shape(1);
  // writer.write_scalar(n_advected_variables, "n_advected_variables");

  // auto advected_labels = numbered_labels("mq%d", n_advected_variables);
  // advected_variables.split_save(writer, advected_labels);
}

// void AllVariables::load(HDF5Reader &reader,
//                         const std::vector<std::string> &labels) {
//   // conserved variables
//   conserved_variables.split_load(reader, labels);

//   // advected variables
//   int n_advected_variables =
//   reader.read_scalar<int>("n_advected_variables");

//   auto advected_labels = numbered_labels("mq%d", n_advected_variables);
//   advected_variables.split_load(reader, advected_labels);
// }

} // namespace zisa
#endif /* end of include guard: ALL_VARIABLES_CPP_WLARDDKM */
