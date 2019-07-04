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

WENOPoly LSQSolver::solve(const array<double, 2, row_major> &rhs) const {

  const auto &x_center = grid->cell_centers(i_cell);
  double length = grid->characteristic_length(i_cell);

  if (order == 1) {
    return WENOPoly{0, {0.0}, x_center, length};
  }

  assert(rhs.size() > 0);

  const auto &moments = grid->normalized_moments(i_cell);
  WENOPoly poly(order - 1, moments, x_center, length);
  constexpr int_t n_vars = WENOPoly::n_vars();

  using RowMajorMatrix
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  Eigen::Index n_coeffs = Eigen::Index(poly.dof(order - 1) - 1);
  Eigen::Map<RowMajorMatrix> coeffs(
      poly.coeffs_ptr() + n_vars, n_coeffs, n_vars);

  coeffs = qr.solve(
      Eigen::Map<const RowMajorMatrix>(rhs.raw(), qr.rows(), n_vars));

  return poly;
}

Eigen::MatrixXd assemble_weno_ao_matrix(const Grid &grid,
                                        const Stencil &stencil) {

  int order = stencil.order();
  double factor = stencil.overfit_factor();

  LOG_ERR_IF(order <= 0,
             string_format("A non-positive convergence order? [%d]", order));

  if (order == 1) {
    return Eigen::MatrixXd::Ones(1, 1);
  }

  auto degree = order - 1;
  auto n_rows = required_stencil_size(degree, factor) - 1;
  auto n_cols = poly_dof<2>(degree) - 1;

  auto A = Eigen::MatrixXd(n_rows, n_cols);

  auto i0 = stencil.global(0);
  auto x0 = grid.cell_centers(i0);
  auto l0 = grid.characteristic_length(i0);
  const auto &C0 = grid.normalized_moments(i0);

  assert(n_rows <= std::numeric_limits<Eigen::Index>::max());
  auto eint = [](zisa::int_t i) {
    assert(i <= std::numeric_limits<Eigen::Index>::max());
    return Eigen::Index(i);
  };

  for (int_t ii = 0; ii < n_rows; ++ii) {
    auto ii_ = eint(ii);
    auto j = stencil.global(ii + 1);
    auto [x_10, x_01, _] = XYZ((grid.cell_centers(j) - x0) / l0);

    auto lj = grid.characteristic_length(j) / l0;
    const auto &Cj = grid.normalized_moments(j);

    if (order >= 2) {
      A(ii_, 0) = x_10;
      A(ii_, 1) = x_01;
    }

    if (order >= 3) {
      auto i_20 = poly_index(2, 0);
      auto i_11 = poly_index(1, 1);
      auto i_02 = poly_index(0, 2);

      auto x_20 = x_10 * x_10;
      auto x_11 = x_10 * x_01;
      auto x_02 = x_01 * x_01;

      auto lj_2 = lj * lj;

      A(ii_, eint(i_20 - 1)) = x_20 - C0(i_20) + lj_2 * Cj(i_20);
      A(ii_, eint(i_11 - 1)) = x_11 - C0(i_11) + lj_2 * Cj(i_11);
      A(ii_, eint(i_02 - 1)) = x_02 - C0(i_02) + lj_2 * Cj(i_02);

      if (order >= 4) {
        auto i_30 = poly_index(3, 0);
        auto i_21 = poly_index(2, 1);
        auto i_12 = poly_index(1, 2);
        auto i_03 = poly_index(0, 3);

        auto x_30 = x_20 * x_10;
        auto x_21 = x_20 * x_01;
        auto x_12 = x_11 * x_01;
        auto x_03 = x_02 * x_01;

        auto lj_3 = lj_2 * lj;

        A(ii_, eint(i_30 - 1))
            = x_30 - C0(i_30) + 3.0 * x_10 * lj_2 * Cj(i_20) + lj_3 * Cj(i_30);

        A(ii_, eint(i_21 - 1)) = x_21 - C0(i_21) + x_01 * lj_2 * Cj(i_20)
                                 + 2.0 * x_10 * lj_2 * Cj(i_11)
                                 + lj_3 * Cj(i_21);

        A(ii_, eint(i_12 - 1)) = x_12 - C0(i_12) + x_10 * lj_2 * Cj(i_02)
                                 + 2.0 * x_01 * lj_2 * Cj(i_11)
                                 + lj_3 * Cj(i_12);

        A(ii_, eint(i_03 - 1))
            = x_03 - C0(i_03) + 3.0 * x_01 * lj_2 * Cj(i_02) + lj_3 * Cj(i_03);

        if (order >= 5) {
          auto i_40 = poly_index(4, 0);
          auto i_31 = poly_index(3, 1);
          auto i_22 = poly_index(2, 2);
          auto i_13 = poly_index(1, 3);
          auto i_04 = poly_index(0, 4);

          auto x_40 = x_30 * x_10;
          auto x_31 = x_30 * x_01;
          auto x_22 = x_21 * x_01;
          auto x_13 = x_12 * x_01;
          auto x_04 = x_03 * x_01;

          auto lj_4 = lj_3 * lj;

          A(ii_, eint(i_40 - 1))
              = x_40 - C0(i_40) + 6.0 * x_20 * lj_2 * Cj(i_20)
                + 4.0 * x_10 * lj_3 * Cj(i_30) + lj_4 * Cj(i_40);

          A(ii_, eint(i_31 - 1))
              = x_31 - C0(i_31) + 3 * x_11 * lj_2 * Cj(i_20)
                + 3.0 * x_20 * lj_2 * Cj(i_11) + x_01 * lj_3 * Cj(i_30)
                + 3.0 * x_10 * lj_3 * Cj(i_21) + lj_4 * Cj(i_31);

          A(ii_, eint(i_22 - 1))
              = x_22 - C0(i_22) + x_02 * lj_2 * Cj(i_20)
                + x_20 * lj_2 * Cj(i_02) + 4 * x_11 * lj_2 * Cj(i_11)
                + 2.0 * x_01 * lj_3 * Cj(i_21) + 2 * x_10 * lj_3 * Cj(i_12)
                + lj_4 * Cj(i_22);

          A(ii_, eint(i_13 - 1))
              = x_13 - C0(i_13) + 3 * x_11 * lj_2 * Cj(i_02)
                + 3.0 * x_02 * lj_2 * Cj(i_11) + x_10 * lj_3 * Cj(i_03)
                + 3.0 * x_01 * lj_3 * Cj(i_12) + lj_4 * Cj(i_13);

          A(ii_, eint(i_04 - 1))
              = x_04 - C0(i_04) + 6.0 * x_02 * lj_2 * Cj(i_02)
                + 4.0 * x_01 * lj_3 * Cj(i_03) + lj_4 * Cj(i_04);
        }
      }
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
