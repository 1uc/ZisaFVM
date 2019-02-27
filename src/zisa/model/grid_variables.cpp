#include <zisa/model/grid_variables.hpp>

namespace zisa {

void save(HDF5Writer &writer,
          const GridVariables &vars,
          const std::vector<std::string> &labels) {
  array<double, 1> component({vars.shape(0)}, device_type::cpu);
  for (int_t k = 0; k < labels.size(); ++k) {
    for (int_t i = 0; i < vars.shape(0); ++i) {
      component[i] = vars(i, k);
    }

    zisa::save(writer, component, labels.at(k));
  }
}

void GridVariables::split_load(HDF5Reader &reader,
                               const std::vector<std::string> &labels) {

  auto dims = reader.dims(labels[0]);
  assert(dims.size() == 1);

  array<double, 1> component(shape_t<1>{dims[0]});

  for (int_t k = 0; k < labels.size(); ++k) {
    zisa::load(reader, component, labels[k]);

    for (int_t i = 0; i < this->shape(1); ++i) {
      (*this)(i, k) = component[i];
    }
  }
}

} // namespace zisa
