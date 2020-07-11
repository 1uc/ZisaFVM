#ifndef GLOBAL_RECONSTRUCTION_IMPL_H_1EKED
#define GLOBAL_RECONSTRUCTION_IMPL_H_1EKED

#include "global_reconstruction_decl.hpp"
#include <zisa/model/all_variables.hpp>
#include <zisa/parallelization/omp.h>
#include <zisa/utils/indent_block.hpp>
#include <zisa/utils/to_string.hpp>
#include <zisa/utils/type_name.hpp>

namespace zisa {

template <class Equilibrium, class RC, class Scaling>
EulerGlobalReconstruction<Equilibrium, RC, Scaling>::EulerGlobalReconstruction(
    const HybridWENOParams &params, array<lrc_t, 1> rc_)
    : params(params),
      rc(std::move(rc_)),
      qbar_allocator(std::make_shared<block_allocator<array<cvars_t, 1>>>(128)),
      tracer_allocator(
          std::make_shared<block_allocator<array<double, 2, column_major>>>(
              128)),
      polys_allocator(
          std::make_shared<block_allocator<array<WENOPoly, 1>>>(128)),
      rhs_allocator(
          std::make_shared<block_allocator<array<double, 2, row_major>>>(128)) {

  n_polys = params.linear_weights.size();
  max_stencil_size = 0;

  for (int_t i = 0; i < rc.size(); ++i) {
    max_stencil_size
        = zisa::max(rc[i].combined_stencil_size(), max_stencil_size);
  }
}

// template <class Equilibrium, class RC, class Scaling>
// EulerGlobalReconstruction<Equilibrium, RC,
// Scaling>::EulerGlobalReconstruction(
//    std::shared_ptr<Grid> grid,
//    const array<StencilFamily, 1> &stencils,
//    const HybridWENOParams &params,
//    const Equilibrium &eq,
//    const Scaling &scaling)
//    : params(params),
//      rc(shape_t<1>{grid->n_cells}),
//      qbar_allocator(std::make_shared<block_allocator<array<cvars_t,
//      1>>>(128)), tracer_allocator(
//          std::make_shared<block_allocator<array<double, 2, column_major>>>(
//              128)),
//      polys_allocator(
//          std::make_shared<block_allocator<array<WENOPoly, 1>>>(128)),
//      rhs_allocator(
//          std::make_shared<block_allocator<array<double, 2, row_major>>>(128))
//          {
//
//  n_polys = params.linear_weights.size();
//  max_stencil_size = 0;
//
//  auto o1_params = HybridWENOParams(
//      {{1}, {"c"}, {1.0}}, {1.0}, params.epsilon, params.exponent);
//
//  for (int_t i = 0; i < grid->n_cells; ++i) {
//    if (stencils[i].size() == 1 && stencils[i].order() == 1) {
//      rc[i] = LocalReconstruction<Equilibrium, RC, Scaling>(
//          grid,
//          LocalEquilibrium(eq),
//          RC(grid, stencils[i], i, o1_params),
//          i,
//          scaling);
//
//    } else {
//      rc[i] = LocalReconstruction<Equilibrium, RC, Scaling>(
//          grid,
//          LocalEquilibrium(eq),
//          RC(grid, stencils[i], i, params),
//          i,
//          scaling);
//    }
//
//    max_stencil_size
//        = zisa::max(rc[i].combined_stencil_size(), max_stencil_size);
//  }
//}

template <class Equilibrium, class RC, class Scaling>
const LocalReconstruction<Equilibrium, RC, Scaling> &
EulerGlobalReconstruction<Equilibrium, RC, Scaling>::operator()(int_t i) const {
  return rc(i);
}

template <class Equilibrium, class RC, class Scaling>
euler_var_t EulerGlobalReconstruction<Equilibrium, RC, Scaling>::operator()(
    int_t i, const XYZ &x) const {
  return rc(i)(x);
}

template <class Equilibrium, class RC, class Scaling>
void EulerGlobalReconstruction<Equilibrium, RC, Scaling>::compute(
    const AllVariables &current_state) {
  auto n_cells = current_state.cvars.shape(0);

#if ZISA_HAS_OPENMP == 1
#pragma omp parallel
#endif
  {
    auto qbar_local = qbar_allocator->allocate(shape_t<1>{max_stencil_size});
    auto tracer_local = tracer_allocator->allocate(
        shape_t<2>{max_stencil_size, current_state.avars.shape(1)});
    auto polys = polys_allocator->allocate(shape_t<1>{n_polys});
    auto rhs = rhs_allocator->allocate(
        shape_t<2>{max_stencil_size, WENOPoly::n_vars()});

#if ZISA_HAS_OPENMP == 1
#pragma omp for ZISA_OMP_FOR_SCHEDULE_DEFAULT
#endif
    for (int_t i = 0; i < n_cells; ++i) {
      set_qbar_local(*qbar_local, current_state, i);
      rc[i].compute(*rhs, *polys, *qbar_local);

      auto tracer_polys = array_view<ScalarPoly, 1>(
          shape_t<1>(polys->shape(0)), (ScalarPoly *)(polys->raw()));
      set_tracer_local(*tracer_local, current_state, i);
      rc[i].compute_tracer(*rhs, tracer_polys, *tracer_local);
    }
  }
}

template <class Equilibrium, class RC, class Scaling>
void EulerGlobalReconstruction<Equilibrium, RC, Scaling>::set_qbar_local(
    array<cvars_t, 1> &qbar_local, const AllVariables &current_state, int_t i) {
  const auto &l2g = rc[i].local2global();
  const auto &cvars = current_state.cvars;

  for (int_t ii = 0; ii < l2g.size(); ++ii) {
    qbar_local(ii) = cvars(l2g[ii]);
  }
}

template <class Equilibrium, class RC, class Scaling>
void EulerGlobalReconstruction<Equilibrium, RC, Scaling>::set_tracer_local(
    array<double, 2, column_major> &tracer_local,
    const AllVariables &current_state,
    int_t i) {

  const auto &l2g = rc[i].local2global();
  const auto &avars = current_state.avars;

  auto n_vars = avars.shape(1);
  for (int_t ii = 0; ii < l2g.size(); ++ii) {
    for (int_t k_vars = 0; k_vars < n_vars; ++k_vars) {
      tracer_local(ii, k_vars) = avars(l2g[ii], k_vars);
    }
  }
}

template <class Equilibrium, class RC, class Scaling>
std::string EulerGlobalReconstruction<Equilibrium, RC, Scaling>::str() const {
  return string_format("EulerGlobalReconstruction<%s>: \n",
                       type_name<RC>().c_str())
         + indent_block(1, zisa::to_string(params));
}

} // namespace zisa

#endif
