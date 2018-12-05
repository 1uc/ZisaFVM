#include <zisa/boundary/no_boundary_condition.hpp>
#include <zisa/model/all_variables.hpp>

namespace zisa {

void NoBoundaryCondition::apply(AllVariables &, double) {
  return; // do nothing
}

std::string NoBoundaryCondition::str() const {
  return "Do-nothing boundary condition.";
}

} // namespace zisa
