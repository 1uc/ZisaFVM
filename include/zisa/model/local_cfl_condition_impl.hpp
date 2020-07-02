#ifndef CFL_CONDITION_IMPL_H_DQRS3
#define CFL_CONDITION_IMPL_H_DQRS3

#include "local_cfl_condition_decl.hpp"

#include <zisa/loops/reduction/min.hpp>
#include <zisa/parallelization/omp.h>

namespace zisa {

template <class Model>
LocalCFL<Model>::LocalCFL(std::shared_ptr<Grid> grid,
                          std::shared_ptr<Model> model,
                          double cfl_number)
    : grid(std::move(grid)), model(std::move(model)), cfl_number(cfl_number) {}

template <class Model>
double LocalCFL<Model>::operator()(const AllVariables &all_variables) {
  const auto &cvars = all_variables.cvars;

  auto f = [this, &cvars](int_t i) {
    auto u = cvars_t(cvars(i));
    auto xvars = model->eos.xvars(u);
    double ev_max = model->max_eigen_value(u, xvars);
    double dx = grid->inradius(i);

    return dx / ev_max;
  };

  return cfl_number * zisa::reduce::min(cell_indices(*grid), f);
}

} // namespace zisa

#endif
