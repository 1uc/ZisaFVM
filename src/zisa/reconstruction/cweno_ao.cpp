#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>

namespace zisa {

auto CWENO_AO::reconstruct(array<WENOPoly, 1> &polys,
                           const array<cvars_t, 1> &qbar) const
    -> decltype(hybridize(polys)) {

  compute_polys(polys, qbar);

  auto k_high = stencils.highest_order_stencil();
  for (int_t k = 0; k < stencils.size(); ++k) {
    if (k_high != k) {
      polys[k_high] -= linear_weights[k] * polys[k];
    }
  }

  polys[k_high] /= linear_weights[k_high];

  return hybridize(polys);
}

} // namespace zisa
