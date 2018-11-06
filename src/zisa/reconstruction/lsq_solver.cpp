#include <zisa/grid/grid.hpp>
#include <zisa/reconstruction/lsq_solver.hpp>

namespace zisa {

Eigen::MatrixXd assemble_weno_ao_matrix(const Grid &grid,
                                        const Stencil &stencil);

LSQSolver::LSQSolver(const std::shared_ptr<Grid> &grid, const Stencil &stencil)
    : grid(grid), i_cell(stencil.global(0)), order(stencil.order()) {

  assert(grid != nullptr);
  assert(stencil.size() > 0);

  auto A = assemble_weno_ao_matrix(*grid, stencil);
  qr.compute(A);
}

auto LSQSolver::solve(const array<double, 1> &rhs) const -> Poly2D<MAX_DEGREE> {

  assert(rhs.size() > 0);

  Eigen::VectorXd coeffs
      = qr.solve(Eigen::Map<const Eigen::VectorXd>(rhs.raw(), qr.rows()));

  const auto &moments = grid->normalized_moments(i_cell);

  if (order == 2) {
    return {{0.0, coeffs(0), coeffs(1)}, {0.0, 0.0, 0.0}};
  }

  if (order == 3) {
    auto i20 = poly_index(2, 0);
    auto i11 = poly_index(1, 1);
    auto i02 = poly_index(0, 2);

    return {{0.0, coeffs(0), coeffs(1), coeffs(2), coeffs(3), coeffs(4)},
            {0.0, 0.0, 0.0, moments(i20), moments(i11), moments(i02)}};
  }

  LOG_ERR("Implement first.");
}

Eigen::MatrixXd assemble_weno_ao_matrix(const Grid &grid,
                                        const Stencil &stencil) {

  int order = stencil.order();
  double factor = stencil.overfit_factor();

  LOG_ERR_IF(order <= 0,
             string_format("A non-positive convergence order? [%d]", order));

  if (order == 0) {
    return Eigen::MatrixXd::Ones(1, 1);
  }

  auto degree = order - 1;
  auto n_rows = required_stencil_size(degree, factor) - 1;
  auto n_cols = poly_dof(degree) - 1;

  auto A = Eigen::MatrixXd(n_rows, n_cols);

  auto i0 = stencil.global(0);
  auto tri0 = grid.triangle(i0);
  auto x0 = grid.cell_centers(i0);
  auto l0 = circum_radius(tri0);
  const auto &C0 = grid.normalized_moments(i0);

  for (int_t ii = 0; ii < n_rows; ++ii) {
    auto j = stencil.global(ii + 1);
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

bool LSQSolver::operator==(const LSQSolver &other) const {
  if (grid != other.grid) {
    return false;
  }

  if (i_cell != other.i_cell) {
    return false;
  }

  if (order != other.order) {
    return false;
  }

  return qr.matrixQR() == other.qr.matrixQR();
}

bool LSQSolver::operator!=(const LSQSolver &other) const {
  return !((*this) == other);
}

} // namespace zisa
