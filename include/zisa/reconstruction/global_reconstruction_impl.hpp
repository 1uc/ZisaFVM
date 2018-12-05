#ifndef GLOBAL_RECONSTRUCTION_IMPL_H_1EKED
#define GLOBAL_RECONSTRUCTION_IMPL_H_1EKED

#include "global_reconstruction_decl.hpp"
#include <zisa/utils/indent_block.hpp>
#include <zisa/utils/to_string.hpp>
#include <zisa/utils/type_name.hpp>

namespace zisa {

template <class RC>
GlobalReconstruction<RC>::GlobalReconstruction(std::shared_ptr<Grid> grid,
                                               const HybridWENO_Params &params,
                                               int_t n_vars)
    : params(params),
      rc(shape_t<1>{grid->n_cells}),
      polys(shape_t<2>{grid->n_cells, n_vars}) {

  int_t max_stencil_size = 0;

  for (int_t i = 0; i < grid->n_cells; ++i) {
    rc[i] = RC(grid, i, params);

    max_stencil_size
        = zisa::max(rc[i].combined_stencil_size(), max_stencil_size);
  }

  qbar_local = array<double, 1>(shape_t<1>{max_stencil_size});
}

template <class RC>
const WENOPoly &GlobalReconstruction<RC>::operator()(int_t i, int_t k) const {
  return polys(i, k);
}

template <class RC>
void GlobalReconstruction<RC>::compute(const AllVariables &current_state) {
  auto n_cells = current_state.cvars.shape(0);
  auto n_vars = current_state.cvars.shape(1);

  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < n_vars; ++k) {
      set_qbar_local(current_state, i, k);
      polys(i, k) = rc[i].reconstruct(qbar_local);
    }
  }
}

template <class RC>
void GlobalReconstruction<RC>::set_qbar_local(const AllVariables &current_state,
                                              int_t i,
                                              int_t k) {
  const auto &l2g = rc[i].local2global();

  for (int_t ii = 0; ii < l2g.size(); ++ii) {
    qbar_local[ii] = current_state.cvars(l2g[ii], k);
  }
}

template <class RC>
std::string GlobalReconstruction<RC>::str() const {
  return string_format("GlobalReconstruction<%s>: \n", type_name<RC>().c_str())
         + indent_block(1, zisa::to_string(params));
}

} // namespace zisa

#endif
