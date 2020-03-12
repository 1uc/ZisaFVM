#include <zisa/boundary/frozen_boundary_condition.hpp>

#include <algorithm>

namespace zisa {

FrozenBC::FrozenBC(const zisa::Grid &grid,
                         const zisa::AllVariables &all_vars) {
  auto n_cells = all_vars.cvars.shape(0);
  auto n_cvars = all_vars.cvars.shape(1);
  auto n_bc_cells = count_ghost_cells(grid);

  indices = array<int_t, 1>(shape_t<1>{n_bc_cells});
  values = array<double, 2>(shape_t<2>{n_bc_cells, n_cvars});

  int_t ii = 0;
  for (int_t i = 0; i < n_cells; ++i) {
    if (grid.cell_flags(i).ghost_cell) {
      indices(ii) = i;
      for (int_t k = 0; k < n_cvars; ++k) {
        values(ii, k) = all_vars.cvars(i, k);
      }
      ++ii;
    }
  }
}

void FrozenBC::apply(AllVariables &u, double) {
  zisa::for_each(flat_range(indices), [this, &u](int_t ii) {
    auto n_cvars = values.shape(1);
    int_t i = indices(ii);

    for (int_t k = 0; k < n_cvars; ++k) {
      u.cvars(i, k) = values(ii, k);
    }
  });
}

std::string FrozenBC::str() const {
  return "Frozen ghost-cell boundary conditions.";
}

int_t FrozenBC::count_ghost_cells(const Grid &grid) const {
  return integer_cast<int_t>(
      std::count_if(grid.cell_flags.begin(),
                    grid.cell_flags.end(),
                    [](const CellFlags &flags) { return flags.ghost_cell; }));
}

}