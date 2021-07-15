// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef CFL_CONDITION_IMPL_H_DQRS3
#define CFL_CONDITION_IMPL_H_DQRS3

#include "local_cfl_condition_decl.hpp"

#include <zisa/loops/reduction/min.hpp>
#include <zisa/parallelization/omp.h>

namespace zisa {

template <class EOS>
LocalCFL<EOS>::LocalCFL(std::shared_ptr<Grid> grid,
                        std::shared_ptr<Euler> euler,
                        std::shared_ptr<LocalEOSState<EOS>> local_eos,
                        double cfl_number)
    : grid(std::move(grid)),
      euler(std::move(euler)),
      local_eos(std::move(local_eos)),
      cfl_number(cfl_number) {}

template <class EOS>
double LocalCFL<EOS>::operator()(const AllVariables &all_variables) {
  const auto &cvars = all_variables.cvars;

  auto f = [this, &cvars](int_t i) {
    const auto &eos = (*local_eos)(i);

    auto u = cvars_t(cvars(i));
    auto xvars = eos->xvars(u);
    double ev_max = euler->max_eigen_value(u, xvars.a);
    double dx = grid->inradius(i);

    return dx / ev_max;
  };

  return cfl_number * zisa::reduce::min(cell_indices(*grid), f);
}

} // namespace zisa

#endif
