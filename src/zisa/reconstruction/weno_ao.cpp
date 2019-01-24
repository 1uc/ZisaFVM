#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/weno_ao.hpp>

namespace zisa {

auto WENO_AO::reconstruct(const array<double, 2> &qbar) const
    -> decltype(hybridize()) {

  compute_polys(qbar);
  return hybridize();
}

} // namespace zisa
