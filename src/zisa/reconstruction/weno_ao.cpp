#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/weno_ao.hpp>

namespace zisa {

auto WENO_AO::reconstruct(array<WENOPoly, 1> &polys,
                          const array<cvars_t, 1> &qbar) const
    -> decltype(hybridize(polys)) {

  compute_polys(polys, qbar);
  return hybridize(polys);
}

} // namespace zisa
