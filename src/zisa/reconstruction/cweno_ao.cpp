#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>

namespace zisa {

CWENO_AO::CWENO_AO(const std::shared_ptr<Grid> &grid,
                   int_t i_cell,
                   const HybridWENO_Params &params)
    : super(grid, i_cell, params) {
  k_high = highest_order_central_stencil(stencils);
}

auto CWENO_AO::reconstruct(const array<double, 1> &qbar) const
    -> decltype(hybridize()) {

  compute_polys(qbar);

  for (int_t k = 0; k < stencils.size(); ++k) {
    if (k_high != k) {
      polys[k_high] -= linear_weights[k] / linear_weights[k_high] * polys[k];
    }
  }

  return hybridize();
}


bool CWENO_AO::operator==(const CWENO_AO &other) const {
  if (k_high != other.k_high) {
    return false;
  }

  return static_cast<const HybridWENO &>(*this)
         == static_cast<const HybridWENO &>(other);
}

bool CWENO_AO::operator!=(const CWENO_AO &other) const {
  return !(*this == other);
}

} // namespace zisa
