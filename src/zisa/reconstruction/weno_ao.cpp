#include <algorithm>
#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/weno_ao.hpp>

namespace zisa {

WENO_AO::WENO_AO(const std::shared_ptr<Grid> &grid, int_t i_cell)
    : n_stencils(4),
      grid(grid),
      i_cell(i_cell),
      stencils(shape_t<1>{n_stencils}),
      qr(shape_t<1>{n_stencils}) {

  compute_stencils();
  if (order >= 2) {
    compute_qr();
  }

  rhs = array<double, 1>(shape_t<1>{required_stencil_size(order - 1, 2.0)});
}

auto WENO_AO::reconstruct(const LocalBuffer &qbar) const
    -> Poly2D<max_degree()> {

  if (order == 1) {
    return {{qbar(0)}, {0.0}};
  }

  int_t k = 0;
  for (int_t i = 0; i < qr(0).rows(); ++i) {
    rhs(i) = qbar(stencils(k)(i + 1).local) - qbar(0);
  }

  Eigen::VectorXd coeffs
      = qr(k).solve(Eigen::Map<const Eigen::VectorXd>(rhs.raw(), qr(k).rows()));

  return convert_to_poly(qbar, coeffs);
}

Poly2D<WENO_AO::max_degree()>
WENO_AO::convert_to_poly(const array<double, 1> &qbar,
                         const Eigen::VectorXd &coeffs) const {
  const auto &moments = grid->normalized_moments(i_cell);
  if (order == 2) {
    return {{qbar(0), coeffs(0), coeffs(1)}, {0.0, 0.0, 0.0}};
  }

  if (order == 3) {
    auto i20 = poly_index(2, 0);
    auto i11 = poly_index(1, 1);
    auto i02 = poly_index(0, 2);

    return {{qbar(0), coeffs(0), coeffs(1), coeffs(2), coeffs(3), coeffs(4)},
            {0.0, 0.0, 0.0, moments(i20), moments(i11), moments(i02)}};
  }

  LOG_ERR("Implement first.");
}

void WENO_AO::compute_qr() {
  for (int_t i = 0; i < stencils.shape(0); ++i) {
    auto A = assemble_weno_ao_matrix(
        *grid, stencils(i), order, (i == 0 ? 2.0 : 1.5));

    qr(i).compute(A);
  }
}

Eigen::MatrixXd assemble_weno_ao_matrix(const Grid &grid,
                                        const array<LocalIndex, 1> &stencil,
                                        int order,
                                        double factor) {
  LOG_ERR_IF(order <= 1,
             string_format("You should not be calling this. [%d]", order));

  auto degree = order - 1;
  auto n_rows = required_stencil_size(degree, factor) - 1;
  auto n_cols = poly_dof(degree) - 1;

  auto A = Eigen::MatrixXd(n_rows, n_cols);

  auto i0 = stencil(0).global;
  auto tri0 = grid.triangle(i0);
  auto x0 = grid.cell_centers(i0);
  auto l0 = circum_radius(tri0);
  const auto &C0 = grid.normalized_moments(i0);

  for (int_t ii = 0; ii < n_rows; ++ii) {
    auto j = stencil(ii + 1).global;
    auto trij = grid.triangle(j);
    XY xj = XY((grid.cell_centers(j) - x0) / l0);

    auto lj = circum_radius(trij) / l0;
    const auto &Cj = grid.normalized_moments(j);

    if (order >= 2) {
      A(ii, 0) = xj(0);
      A(ii, 1) = xj(1);
    }

    if (order >= 3) {
      auto i_20 = poly_index(2, 0);
      auto i_11 = poly_index(1, 1);
      auto i_02 = poly_index(0, 2);

      A(ii, i_20 - 1) = xj(0) * xj(0) - C0(i_20) + lj * lj * Cj(i_20);
      A(ii, i_11 - 1) = xj(0) * xj(1) - C0(i_11) + lj * lj * Cj(i_11);
      A(ii, i_02 - 1) = xj(1) * xj(1) - C0(i_02) + lj * lj * Cj(i_02);
    }
    if (order == 4) {
      LOG_ERR("Implement first.");
    }

    if (order == 5) {
      LOG_ERR("Derive & implement first.");
    }
    if (order >= 6) {
      LOG_ERR("Good luck. :) ");
    }
  }

  return A;
}

void WENO_AO::compute_stencils() {
  int_t max_points = required_stencil_size(max_degree(), 2.0);

  auto c = zisa::central_stencil(*grid, i_cell, max_points);

  auto b0 = zisa::biased_stencil(*grid, i_cell, 0, max_points);
  auto b1 = zisa::biased_stencil(*grid, i_cell, 1, max_points);
  auto b2 = zisa::biased_stencil(*grid, i_cell, 2, max_points);

  order = zisa::min(deduce_max_order(c, 2.0),
                    deduce_max_order(b0, 1.5),
                    deduce_max_order(b1, 1.5),
                    deduce_max_order(b2, 1.5));

  order = zisa::min(max_degree() + 1, order);

  c.resize(required_stencil_size(order - 1, 2.0));
  b0.resize(required_stencil_size(order - 1, 1.5));
  b1.resize(required_stencil_size(order - 1, 1.5));
  b2.resize(required_stencil_size(order - 1, 1.5));

  l2g.reserve(c.size());

  stencils(0) = assign_local_indices(c, l2g);
  stencils(1) = assign_local_indices(b0, l2g);
  stencils(2) = assign_local_indices(b1, l2g);
  stencils(3) = assign_local_indices(b2, l2g);
}

array<LocalIndex, 1>
WENO_AO::assign_local_indices(const std::vector<int_t> &global_indices,
                              std::vector<int_t> &l2g) {

  assert(global_indices.size() > 0);

  auto a = array<LocalIndex, 1>(shape_t<1>{global_indices.size()});

  for (int_t i = 0; i < global_indices.size(); ++i) {
    a(i).global = global_indices[i];

    auto local_index_ptr = std::find(l2g.begin(), l2g.end(), a(i).global);
    a(i).local = local_index_ptr - l2g.begin();

    if (local_index_ptr == l2g.end()) {
      l2g.push_back(a(i).global);
    }
  }

  return a;
}

const std::vector<int_t> &WENO_AO::local2global() const { return l2g; }

int deduce_max_order(const std::vector<int_t> &stencil, double factor) {
  auto stencil_size = stencil.size();

  int deg = 0;
  while (required_stencil_size(deg + 1, factor) <= stencil_size) {
    deg += 1;
  }

  return deg + 1; // order is the degree plus one
}

int_t required_stencil_size(int deg, double factor) {
  if (deg == 0) {
    return 1;
  }

  return int_t(double(poly_dof(deg) - 1) * factor);
}

std::vector<int_t>
central_stencil(const Grid &grid, int_t i_center, int_t n_points) {

  return biased_stencil(
      grid, i_center, n_points, Cone({0.0, 0.0}, {1.0, 0.0}, -1.1));
}
std::vector<int_t>
biased_stencil(const Grid &grid, int_t i_center, int_t k, int_t n_points) {

  const auto &x_center = grid.cell_centers(i_center);
  const auto &v0 = grid.vertex(i_center, k);
  const auto &v1 = grid.vertex(i_center, (k + 1) % grid.max_neighbours);

  return biased_stencil(grid, i_center, n_points, Cone(v0, x_center, v1));
}

std::vector<int_t> biased_stencil(const Grid &grid,
                                  int_t i_center,
                                  int_t n_points,
                                  const Cone &cone) {

  int_t max_points = 5 * n_points;
  int_t max_neighbours = grid.max_neighbours;

  std::vector<int_t> candidates;
  candidates.reserve(3 * max_points);
  candidates.push_back(i_center);

  auto not_found = [&candidates](int_t cand) {
    return std::find(candidates.begin(), candidates.end(), cand)
           == candidates.end();
  };

  auto is_inside = [&cone, &grid](int_t cand) {
    return cone.is_inside(grid.cell_centers(cand));
  };

  for (int_t p = 0; p < max_points; ++p) {
    if (p >= candidates.size()) {
      break;
    }

    int_t j = candidates[p];
    for (int_t k = 0; k < max_neighbours; ++k) {
      if (!grid.is_valid(j, k)) {
        continue;
      }

      int_t cand = grid.neighbours(j, k);

      if (not_found(cand) && is_inside(cand)) {
        candidates.push_back(cand);
      }
    }
  }

  auto closer = [&grid, i_center](int_t i1, int_t i2) {
    const auto &x_center = grid.cell_centers(i_center);
    return zisa::norm(grid.cell_centers(i1) - x_center)
           < zisa::norm(grid.cell_centers(i2) - x_center);
  };

  std::sort(candidates.begin(), candidates.end(), closer);
  candidates.resize(zisa::min(n_points, candidates.size()));
  return candidates;
}

} // namespace zisa
