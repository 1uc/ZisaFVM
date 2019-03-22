#ifndef GLOBAL_RECONSTRUCTION_IMPL_H_1EKED
#define GLOBAL_RECONSTRUCTION_IMPL_H_1EKED

#include "global_reconstruction_decl.hpp"
#include <zisa/model/all_variables.hpp>
#include <zisa/parallelization/omp.h>
#include <zisa/utils/indent_block.hpp>
#include <zisa/utils/to_string.hpp>
#include <zisa/utils/type_name.hpp>

namespace zisa {

template <class Equilibrium, class RC>
EulerGlobalReconstruction<Equilibrium, RC>::EulerGlobalReconstruction(
    std::shared_ptr<Grid> grid,
    const HybridWENOParams &params,
    const Equilibrium &eq)
    : params(params),
      rc(shape_t<1>{grid->n_cells}),
      allocator(std::make_unique<block_allocator<array<cvars_t, 1>>>(128)) {

  max_stencil_size = 0;

  for (int_t i = 0; i < grid->n_cells; ++i) {
    rc[i] = LocalReconstruction<Equilibrium, RC>(
        grid, {eq}, {grid, i, params}, grid->triangle(i));

    max_stencil_size
        = zisa::max(rc[i].combined_stencil_size(), max_stencil_size);
  }
}

template <class Equilibrium, class RC>
const LocalReconstruction<Equilibrium, RC> &
EulerGlobalReconstruction<Equilibrium, RC>::operator()(int_t i) const {
  return rc(i);
}

template <class Equilibrium, class RC>
euler_var_t EulerGlobalReconstruction<Equilibrium, RC>::
operator()(int_t i, const XYZ &x) const {
  return rc(i)(x);
}

template <class Equilibrium, class RC>
void EulerGlobalReconstruction<Equilibrium, RC>::compute(
    const AllVariables &current_state) {
  auto n_cells = current_state.cvars.shape(0);

#pragma omp parallel
  {
    auto qbar_local = allocator->allocate(shape_t<1>{max_stencil_size});

#pragma omp for ZISA_OMP_FOR_SCHEDULE_DEFAULT
    for (int_t i = 0; i < n_cells; ++i) {
      set_qbar_local(*qbar_local, current_state, i);
      rc[i].compute(*qbar_local);
    }
  }
}

template <class Equilibrium, class RC>
void EulerGlobalReconstruction<Equilibrium, RC>::set_qbar_local(
    array<cvars_t, 1> &qbar_local, const AllVariables &current_state, int_t i) {
  const auto &l2g = rc[i].local2global();
  const auto &cvars = current_state.cvars;

  for (int_t ii = 0; ii < l2g.size(); ++ii) {
    qbar_local(ii) = cvars(l2g[ii]);
  }
}

template <class Equilibrium, class RC>
std::string EulerGlobalReconstruction<Equilibrium, RC>::str() const {
  return string_format("EulerGlobalReconstruction<%s>: \n",
                       type_name<RC>().c_str())
         + indent_block(1, zisa::to_string(params));
}

} // namespace zisa

#endif
