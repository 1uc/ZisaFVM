#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/weno_ao.hpp>

namespace zisa {

WENO_AO::WENO_AO(const std::shared_ptr<Grid> &grid,
                 int_t i_cell,
                 const HybridWENO_Params &params)
    : super(grid, i_cell, params) {}

auto WENO_AO::reconstruct(const array<double, 1> &qbar) const
    -> decltype(hybridize()) {

  compute_polys(qbar);
  return hybridize();
}

} // namespace zisa
