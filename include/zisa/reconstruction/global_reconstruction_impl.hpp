#ifndef GLOBAL_RECONSTRUCTION_IMPL_H_1EKED
#define GLOBAL_RECONSTRUCTION_IMPL_H_1EKED

#include "global_reconstruction_decl.hpp"
#include <zisa/model/all_variables.hpp>
#include <zisa/utils/indent_block.hpp>
#include <zisa/utils/to_string.hpp>
#include <zisa/utils/type_name.hpp>

#include <omp.h>

namespace zisa {

template <class RC>
GlobalReconstruction<RC>::GlobalReconstruction(std::shared_ptr<Grid> grid,
                                               const HybridWENO_Params &params,
                                               int_t n_vars)
    : params(params),
      rc(shape_t<1>{grid->n_cells}),
      polys(shape_t<1>{grid->n_cells}) {

  assert(n_vars == 5);

  int_t max_stencil_size = 0;

  for (int_t i = 0; i < grid->n_cells; ++i) {
    rc[i] = RC(grid, i, params);

    max_stencil_size
        = zisa::max(rc[i].combined_stencil_size(), max_stencil_size);
  }

  int_t n_threads = int_t(omp_get_max_threads());
  for (int_t i = 0; i < n_threads; ++i) {
    qbar_local.push_back(
        array<double, 2>(shape_t<2>{max_stencil_size, n_vars}));
  }
}

template <class RC>
const WENOPoly &GlobalReconstruction<RC>::operator()(int_t i) const {
  return polys(i);
}

template <class RC>
void GlobalReconstruction<RC>::compute(const AllVariables &current_state) {
  auto n_cells = current_state.cvars.shape(0);

#pragma omp parallel for schedule(guided)
  for (int_t i = 0; i < n_cells; ++i) {
    auto thread_id = int_t(omp_get_thread_num());
    set_qbar_local(current_state, i);
    polys(i) = rc[i].reconstruct(qbar_local[thread_id]);
  }
}

template <class RC>
void GlobalReconstruction<RC>::set_qbar_local(const AllVariables &current_state,
                                              int_t i) {
  const auto &l2g = rc[i].local2global();
  auto thread_id = int_t(omp_get_thread_num());
  auto &qbar_local = this->qbar_local[thread_id];

  for (int_t ii = 0; ii < l2g.size(); ++ii) {
    for (int_t k = 0; k < qbar_local.shape(1); ++k) {
      qbar_local(ii, k) = current_state.cvars(l2g[ii], k);
    }
  }
}

template <class RC>
std::string GlobalReconstruction<RC>::str() const {
  return string_format("GlobalReconstruction<%s>: \n", type_name<RC>().c_str())
         + indent_block(1, zisa::to_string(params));
}

} // namespace zisa

#endif
