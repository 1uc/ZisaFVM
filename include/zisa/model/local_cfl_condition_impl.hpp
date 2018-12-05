#ifndef CFL_CONDITION_IMPL_H_DQRS3
#define CFL_CONDITION_IMPL_H_DQRS3

#include "local_cfl_condition_decl.hpp"

namespace zisa {

template <class Model>
LocalCFL<Model>::LocalCFL(std::shared_ptr<Grid> grid,
                          Model model,
                          double cfl_number)
    : grid(grid), model(std::move(model)), cfl_number(cfl_number) {}

template <class Model>
double LocalCFL<Model>::operator()(const AllVariables &all_variables) {
  const auto &u = all_variables.cvars;

  double dt_inv = 0.0;

  for (auto &&[i, tri] : triangles(*grid)) {

    double ev_max = model.max_eigen_value(cvars_t(u(i)));
    double dx = inradius(tri);

    dt_inv = zisa::max(dt_inv, ev_max / dx);
  }

  return cfl_number / dt_inv;
}

} // namespace zisa

#endif
