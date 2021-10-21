// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/grid/grid.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/reconstruction/assemble_weno_ao_matrix.hpp>
#include <zisa/reconstruction/lsq_solver.hpp>

namespace zisa {
LSQSolver::LSQSolver(const std::shared_ptr<Grid> &grid, const Stencil &stencil)
    : grid(grid), i_cell(stencil.global(0)), order(stencil.order()) {

  assert(grid != nullptr);
  assert(stencil.size() > 0);

  A = assemble_weno_ao_matrix(*grid, stencil);
  ldlt.compute(A.transpose() * A);
}

int LSQSolver::n_dims() const { return grid->n_dims(); }

template <class Poly>
Poly LSQSolver::solve_impl(
    const array_const_view<double, 2, row_major> &rhs) const {

  const auto &x_center = grid->cell_centers(i_cell);
  double length = grid->characteristic_length(i_cell);

  int n_dims = grid->n_dims();

  if (order == 1) {
    return Poly(0, {0.0}, x_center, length, n_dims);
  }

  assert(rhs.size() > 0);
  assert(rhs.shape(1) == Poly::n_vars());

  const auto &moments = grid->normalized_moments(i_cell);
  Poly poly(order - 1, moments, x_center, length, n_dims);
  constexpr int_t n_vars = Poly::n_vars();

  using RowMajorMatrix
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  auto n_coeffs = Eigen::Index(poly.dof() - 1);
  auto coeffs = Eigen::Map<RowMajorMatrix>(
      poly.coeffs_ptr() + n_vars, n_coeffs, n_vars);

  auto mapped_rhs
      = Eigen::Map<const RowMajorMatrix>(rhs.raw(), A.rows(), n_vars);

  coeffs = ldlt.solve(A.transpose() * mapped_rhs);

  return poly;
}

template <class Poly>
Poly LSQSolver::solve(const array_const_view<double, 2, row_major> &rhs) const {
  return solve_impl<Poly>(rhs);
}

template WENOPoly LSQSolver::solve<WENOPoly>(
    const array_const_view<double, 2, row_major> &rhs) const;

template ScalarPoly LSQSolver::solve<ScalarPoly>(
    const array_const_view<double, 2, row_major> &rhs) const;

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

  // FIXME
  //  return qr.matrixQR() == other.qr.matrixQR();
  return true;
}

bool LSQSolver::operator!=(const LSQSolver &other) const {
  return !((*this) == other);
}

} // namespace zisa
