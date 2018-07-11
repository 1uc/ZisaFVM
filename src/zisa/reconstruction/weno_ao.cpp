#include <algorithm>
#include <zisa/reconstruction/weno_ao.hpp>

namespace zisa {

WENO_AO::WENO_AO(const Grid &grid, int_t i_cell) : stencils(shape_t<1>{4ul}) {

  // int max_order = 3;
  // int_t n_points = 6 * 2;

  // stencils(0) = central_stencil(grid, i_cell, n_points);
  // stencils(1) = biased_stencil(grid, i_cell, n_points, );
  // stencils(2) = biased_stencil(grid, i_cell, n_points);
  // stencils(3) = biased_stencil(grid, i_cell, n_points);
}

std::vector<int_t> central_stencil(const Grid &grid, int_t i, int_t n_points) {

  int_t max_points = 5 * n_points;
  int_t max_neighbours = grid.max_neighbours;

  std::vector<int_t> candidates;
  candidates.reserve(3 * max_points);
  candidates.push_back(i);

  auto not_found = [&candidates](int_t cand) {
    return std::find(candidates.begin(), candidates.end(), cand)
           == candidates.end();
  };

  for (int_t p = 0; p < max_points; ++p) {
    for (int_t k = 0; k < max_neighbours; ++k) {

      int_t j = candidates[p];
      if (grid.is_valid(j, k)) {
        int_t cand = grid.neighbours(j, k);
        if (not_found(cand)) {
          candidates.push_back(cand);
        }
      }
    }
  }

  auto closer = [&grid, i](int_t i1, int_t i2) {
    return zisa::norm(grid.cell_centers(i1) - grid.cell_centers(i))
           < zisa::norm(grid.cell_centers(i2) - grid.cell_centers(i));
  };

  std::sort(candidates.begin(), candidates.end(), closer);
  candidates.resize(n_points);
  return candidates;
}

Cone::Cone(const XY &x, const XY &dir, double cos_angle)
    : x(x), dir(normalize(dir)), cos_angle(cos_angle) {}

bool Cone::is_inside(const XY &y) const {
  return zisa::dot(dir, y - x) > cos_angle;
}

std::vector<int_t> biased_stencil(
    const Grid &grid, int_t i, int_t n_points, XY dir, double cos_angle) {

  for (int_t k = 0; k < 3; ++k) {
    int_t j = grid.neighbours(i, k);

    const auto &x = grid.cell_centers(j);
  }

  LOG_ERR("Implement first.");
}

} // namespace zisa
