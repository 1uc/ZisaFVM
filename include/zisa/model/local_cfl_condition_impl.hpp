#ifndef CFL_CONDITION_IMPL_H_DQRS3
#define CFL_CONDITION_IMPL_H_DQRS3

#include "local_cfl_condition_decl.hpp"

#include <zisa/parallelization/omp.h>

namespace zisa {

template <class Model>
LocalCFL<Model>::LocalCFL(std::shared_ptr<Grid> grid,
                          std::shared_ptr<Model> model,
                          double cfl_number)
    : grid(grid), model(std::move(model)), cfl_number(cfl_number) {}

template <class Model>
double LocalCFL<Model>::operator()(const AllVariables &all_variables) {
  const auto &u = all_variables.cvars;

  double dt_inv = 0.0;

  int_t n_cells = grid->n_cells;

#pragma omp parallel for ZISA_OMP_FOR_SCHEDULE_DEFAULT
  for (int_t i = 0; i < n_cells; ++i) {
    auto tri = grid->triangle(i);

    double ev_max = model->max_eigen_value(cvars_t(u(i)));
    double dx = inradius(tri);

    dt_inv = zisa::max(dt_inv, ev_max / dx);
  }

  return cfl_number / dt_inv;
}

} // namespace zisa

#endif
