#include <zisa/model/grid_variables.hpp>

namespace zisa {

void save(HierarchicalWriter &writer,
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

GridVariables GridVariables::load(HierarchicalReader &reader,
                                  const std::vector<std::string> &labels) {

  if (labels.empty()) {
    return GridVariables{};
  }

  auto n_cells = int_t(reader.dims(labels[0])[0]);
  auto n_vars = int_t(labels.size());
  auto shape = shape_t<2>{n_cells, n_vars};

  GridVariables vars(shape, device_type::cpu);
  GridVariables::load(reader, vars, labels);
  return vars;
}

void GridVariables::load(HierarchicalReader &reader,
                         GridVariables &vars,
                         const std::vector<std::string> &labels) {

  if (labels.empty()) {
    return;
  }

  auto n_cells = int_t(reader.dims(labels[0])[0]);
  auto n_vars = int_t(labels.size());
  auto shape = shape_t<2>{n_cells, n_vars};

  LOG_ERR_IF(vars.shape(0) < shape(0), "Invalid shape.");
  LOG_ERR_IF(vars.shape(1) != shape(1), "Invalid shape.");

  for (int_t k = 0; k < labels.size(); ++k) {
    auto component = array<double, 1>::load(reader, labels[k]);
    LOG_ERR_IF(component.shape(0) != n_cells, "Reloading non-uniform arrays.");

    for (int_t i = 0; i < component.shape(0); ++i) {
      vars(i, k) = component[i];
    }
  }
}

DereferenceConstGridVariables::DereferenceConstGridVariables(
    const GridVariables &grid_vars)
    : grid_vars(grid_vars) {}

} // namespace zisa
