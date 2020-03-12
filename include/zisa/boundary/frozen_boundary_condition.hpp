#ifndef ZISA_FROZEN_BOUNDARY_CONDITION_HPP
#define ZISA_FROZEN_BOUNDARY_CONDITION_HPP

#include <zisa/boundary/boundary_condition.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/model/all_variables.hpp>

namespace zisa {

class FrozenBC : public BoundaryCondition {
public:
  FrozenBC(const Grid &grid, const AllVariables &all_vars);

  /// Apply the boundary conditions to `u`.
  void apply(AllVariables &u, double /* t */) override;

  std::string str() const override;

protected:
  int_t count_ghost_cells(const Grid &grid) const;

private:
  array<int_t, 1> indices;
  array<double, 2> values;
};
}

#endif // ZISA_FROZEN_BOUNDARY_CONDITION_HPP
