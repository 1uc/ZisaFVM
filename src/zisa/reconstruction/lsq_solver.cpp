#include <zisa/grid/grid.hpp>
#include <zisa/reconstruction/lsq_solver.hpp>

namespace zisa {
Eigen::MatrixXd allocate_weno_ao_matrix(const Grid &grid,
                                        const Stencil &stencil);

Eigen::MatrixXd assemble_weno_ao_matrix(const Grid &grid,
                                        const Stencil &stencil);

void assemble_weno_ao_matrix(Eigen::MatrixXd &A,
                             const Grid &grid,
                             const Stencil &stencil);

void assemble_2d_weno_ao_matrix(Eigen::MatrixXd &A,
                                const Grid &grid,
                                const Stencil &stencil);

void assemble_3d_weno_ao_matrix(Eigen::MatrixXd &A,
                                const Grid &grid,
                                const Stencil &stencil);

LSQSolver::LSQSolver(const std::shared_ptr<Grid> &grid, const Stencil &stencil)
    : grid(grid), i_cell(stencil.global(0)), order(stencil.order()) {

  assert(grid != nullptr);
  assert(stencil.size() > 0);

  auto A = assemble_weno_ao_matrix(*grid, stencil);
  qr.compute(A);
}

int LSQSolver::n_dims() const { return grid->n_dims(); }

WENOPoly LSQSolver::solve(const array<double, 2, row_major> &rhs) const {

  const auto &x_center = grid->cell_centers(i_cell);
  double length = grid->characteristic_length(i_cell);

  int n_dims = grid->n_dims();

  if (order == 1) {
    return WENOPoly(0, {0.0}, x_center, length, n_dims);
  }

  assert(rhs.size() > 0);

  const auto &moments = grid->normalized_moments(i_cell);
  WENOPoly poly(order - 1, moments, x_center, length, n_dims);
  constexpr int_t n_vars = WENOPoly::n_vars();

  using RowMajorMatrix
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  auto n_coeffs = Eigen::Index(poly.dof() - 1);
  auto coeffs = Eigen::Map<RowMajorMatrix>(
      poly.coeffs_ptr() + n_vars, n_coeffs, n_vars);

  auto mapped_rhs
      = Eigen::Map<const RowMajorMatrix>(rhs.raw(), qr.rows(), n_vars);

  coeffs = qr.solve(mapped_rhs);

  return poly;
}

Eigen::MatrixXd assemble_weno_ao_matrix(const Grid &grid,
                                        const Stencil &stencil) {

  auto A = allocate_weno_ao_matrix(grid, stencil);
  assemble_weno_ao_matrix(A, grid, stencil);

  return A;
}

void assemble_weno_ao_matrix(Eigen::MatrixXd &A,
                             const Grid &grid,
                             const Stencil &stencil) {
  int n_dims = grid.n_dims();

  if (n_dims == 2) {
    assemble_2d_weno_ao_matrix(A, grid, stencil);
  } else if (n_dims == 3) {
    assemble_2d_weno_ao_matrix(A, grid, stencil);
    assemble_3d_weno_ao_matrix(A, grid, stencil);
  } else {
    LOG_ERR(string_format("Only works in 2D or 3D. [%d]", n_dims));
  }
}

Eigen::MatrixXd allocate_weno_ao_matrix(const Grid &grid,
                                        const Stencil &stencil) {
  int order = stencil.order();
  int n_dims = grid.n_dims();
  double factor = stencil.overfit_factor();

  LOG_ERR_IF(order <= 0,
             string_format("A non-positive convergence order? [%d]", order));

  if (order == 1) {
    return Eigen::MatrixXd::Ones(1, 1);
  }

  auto degree = order - 1;
  auto n_rows = required_stencil_size(degree, factor, n_dims) - 1;
  auto n_cols = poly_dof(degree, n_dims) - 1;

  return Eigen::MatrixXd(n_rows, n_cols);
}

void assemble_2d_weno_ao_matrix(Eigen::MatrixXd &A,
                                const Grid &grid,
                                const Stencil &stencil) {

  int order = stencil.order();
  if (order == 1) {
    return;
  }

  auto n_dims = grid.n_dims();
  auto n_rows = int_t(A.rows());

  auto i0 = stencil.global(0);
  auto x0 = grid.cell_centers(i0);
  auto l0 = grid.characteristic_length(i0);
  const auto &C0 = grid.normalized_moments(i0);

  assert(n_rows <= std::numeric_limits<Eigen::Index>::max());
  auto eint = [](zisa::int_t i) {
    assert(i <= std::numeric_limits<Eigen::Index>::max());
    return Eigen::Index(i);
  };

  auto idx = [n_dims](int i, int j) {
    return (n_dims == 2) ? poly_index(i, j) : poly_index(i, j, 0);
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
      auto i_20 = idx(2, 0);
      auto i_11 = idx(1, 1);
      auto i_02 = idx(0, 2);

      auto x_20 = x_10 * x_10;
      auto x_11 = x_10 * x_01;
      auto x_02 = x_01 * x_01;

      auto lj_2 = lj * lj;

      A(ii_, eint(i_20 - 1)) = x_20 - C0(i_20) + lj_2 * Cj(i_20);
      A(ii_, eint(i_11 - 1)) = x_11 - C0(i_11) + lj_2 * Cj(i_11);
      A(ii_, eint(i_02 - 1)) = x_02 - C0(i_02) + lj_2 * Cj(i_02);

      if (order >= 4) {
        auto i_30 = idx(3, 0);
        auto i_21 = idx(2, 1);
        auto i_12 = idx(1, 2);
        auto i_03 = idx(0, 3);

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
          auto i_40 = idx(4, 0);
          auto i_31 = idx(3, 1);
          auto i_22 = idx(2, 2);
          auto i_13 = idx(1, 3);
          auto i_04 = idx(0, 4);

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
}

void assemble_3d_weno_ao_matrix(Eigen::MatrixXd &A,
                                const Grid &grid,
                                const Stencil &stencil) {

  int order = stencil.order();
  LOG_ERR_IF(order == 0, "Invalid stencil");
  if (order == 1) {
    return;
  }

  auto n_rows = int_t(A.rows());

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
    auto [x, y, z] = XYZ((grid.cell_centers(j) - x0) / l0);

    auto lj = grid.characteristic_length(j) / l0;
    const auto &Cj = grid.normalized_moments(j);

    if (order >= 2) {
      A(ii_, 2) = z;
    }

    if (order >= 3) {
      auto i_002 = poly_index(0, 0, 2);
      auto i_101 = poly_index(1, 0, 1);
      auto i_011 = poly_index(0, 1, 1);

      auto x_002 = z * z;
      auto x_101 = x * z;
      auto x_011 = y * z;

      auto lj_2 = lj * lj;

      A(ii_, eint(i_002 - 1)) = x_002 - C0(i_002) + lj_2 * Cj(i_002);
      A(ii_, eint(i_101 - 1)) = x_101 - C0(i_101) + lj_2 * Cj(i_101);
      A(ii_, eint(i_011 - 1)) = x_011 - C0(i_011) + lj_2 * Cj(i_011);

      if (order >= 4) {
        auto i_003 = poly_index(0, 0, 3);
        auto i_102 = poly_index(1, 0, 2);
        auto i_012 = poly_index(0, 1, 2);
        auto i_201 = poly_index(2, 0, 1);
        auto i_111 = poly_index(1, 1, 1);
        auto i_021 = poly_index(0, 2, 1);
        auto i_200 = poly_index(2, 0, 0);
        auto i_020 = poly_index(0, 2, 0);
        auto i_110 = poly_index(1, 1, 0);

        auto x_003 = z * z * z;
        auto x_102 = x * z * z;
        auto x_012 = y * z * z;
        auto x_201 = x * x * z;
        auto x_111 = x * y * z;
        auto x_021 = y * y * z;

        auto lj_3 = lj * lj * lj;

        // clang-format off
        A(ii_, eint(i_003 - 1))
            = x_003 - C0(i_003)
            + 3.0 * z * lj_2 * Cj(i_002)
            + lj_3 * Cj(i_003);

        A(ii_, eint(i_102 - 1))
            = x_102 - C0(i_102)
            + x * lj_2 * Cj(i_002)
            + 2.0 * z * lj_2 * Cj(i_101)
            + lj_3 * Cj(i_102);

        A(ii_, eint(i_012 - 1))
            = x_012 - C0(i_012)
            + y * lj_2 * Cj(i_002)
            + 2.0 * z * lj_2 * Cj(i_011)
            + lj_3 * Cj(i_012);

        A(ii_, eint(i_201 - 1))
            = x_201 - C0(i_201)
            + z * lj_2 * Cj(i_200)
            + 2.0 * x * lj_2 * Cj(i_101)
            + lj_3 * Cj(i_201);

        A(ii_, eint(i_021 - 1))
                = x_021 - C0(i_021)
                  + z * lj_2 * Cj(i_020)
                  + 2.0 * y * lj_2 * Cj(i_011)
                  + lj_3 * Cj(i_021);

        A(ii_, eint(i_111 - 1))
            = x_111 - C0(i_111)
            + z * lj_2 * Cj(i_110)
            + y * lj_2 * Cj(i_101)
            + x * lj_2 * Cj(i_011)
            + lj_3 * Cj(i_111);
        // clang-format on
      }
    }
  }
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
