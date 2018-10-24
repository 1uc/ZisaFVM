#include <algorithm>
#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/weno_ao.hpp>

namespace zisa {

WENO_AO::WENO_AO(const std::shared_ptr<Grid> &grid,
                 int_t i_cell,
                 const WENO_AO_Params &params)
    : stencils(grid, i_cell, params.stencil_family_params),
      lsq_solvers(grid, stencils),
      linear_weights(params.linear_weights) {

  rhs = array<double, 1>(shape_t<1>{stencils.combined_stencil_size()});
}

auto WENO_AO::reconstruct(const array<double, 1> &qbar) const
    -> Poly2D<MAX_DEGREE> {

  assert(stencils.size() == 1);

  auto p_avg = Poly2D<MAX_DEGREE>{{qbar(0)}, {0.0}};
  if (stencils.order() == 1) {
    return p_avg;
  }

  int_t k = 0;
  for (int_t i = 0; i < stencils[0].size() - 1; ++i) {
    rhs(i) = qbar(stencils[k].local(i + 1)) - qbar(0);
  }

  return p_avg + lsq_solvers[0].solve(rhs);
}

auto WENO_AO::local2global() const -> decltype(stencils.local2global()) {
  return stencils.local2global();
}

// Poly2D<WENO_AO::max_degree()>
// WENO_AO::convert_to_poly(const array<double, 1> &qbar,
//                          const Eigen::VectorXd &coeffs) const {
//   const auto &moments = grid->normalized_moments(i_cell);
//   if (order == 2) {
//     return {{qbar(0), coeffs(0), coeffs(1)}, {0.0, 0.0, 0.0}};
//   }

//   if (order == 3) {
//     auto i20 = poly_index(2, 0);
//     auto i11 = poly_index(1, 1);
//     auto i02 = poly_index(0, 2);

//     return {{qbar(0), coeffs(0), coeffs(1), coeffs(2), coeffs(3), coeffs(4)},
//             {0.0, 0.0, 0.0, moments(i20), moments(i11), moments(i02)}};
//   }

//   LOG_ERR("Implement first.");
// }

// void WENO_AO::compute_qr() {
//   for (int_t i = 0; i < stencils.shape(0); ++i) {
//     auto A = assemble_weno_ao_matrix(
//         *grid, stencils(i), order, (i == 0 ? 2.0 : 1.5));

//     qr(i).compute(A);
//   }
// }

// void WENO_AO::compute_stencils() {
// LOG_ERR("Implement");
// int_t max_points = required_stencil_size(max_degree(), 2.0);

// auto c = zisa::central_stencil(*grid, i_cell, max_points);

// auto b0 = zisa::biased_stencil(*grid, i_cell, 0, max_points);
// auto b1 = zisa::biased_stencil(*grid, i_cell, 1, max_points);
// auto b2 = zisa::biased_stencil(*grid, i_cell, 2, max_points);

// order = zisa::min(deduce_max_order(c, 2.0),
//                   deduce_max_order(b0, 1.5),
//                   deduce_max_order(b1, 1.5),
//                   deduce_max_order(b2, 1.5));

// order = zisa::min(max_degree() + 1, order);

// c.resize(required_stencil_size(order - 1, 2.0));
// b0.resize(required_stencil_size(order - 1, 1.5));
// b1.resize(required_stencil_size(order - 1, 1.5));
// b2.resize(required_stencil_size(order - 1, 1.5));

// l2g.reserve(c.size());

// stencils(0) = assign_local_indices(c, l2g);
// stencils(1) = assign_local_indices(b0, l2g);
// stencils(2) = assign_local_indices(b1, l2g);
// stencils(3) = assign_local_indices(b2, l2g);
// }

} // namespace zisa
