// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef LSQ_SOLVER_H_7JQG4
#define LSQ_SOLVER_H_7JQG4

#include <Eigen/Dense>

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/stencil.hpp>
#include <zisa/reconstruction/weno_poly.hpp>

namespace zisa {

/// Solve the least-squares problem for the reconstruction.
/** WENO-AO defines the reconstruction polynomial as the polynomial
 * which approximates the cell-averages on a given stencil the best, in a
 * least-squares sense.
 *
 * This class solves those LSQ problems.
 */
class LSQSolver {
private:
  using LDLT = Eigen::LDLT<Eigen::MatrixXd>;

public:
  LSQSolver() = default;
  LSQSolver(const std::shared_ptr<Grid> &grid, const Stencil &stencil);

  /// Solve the LSQ problem with right-hand side `rhs`.
  template <class Poly>
  Poly solve(const array_const_view<double, 2, row_major> &rhs) const;

  /// Indistinguishable by calls to the public interface.
  bool operator==(const LSQSolver &other) const;

  /// Can be distinguished by calls to the public interface.
  bool operator!=(const LSQSolver &other) const;

private:
  template <class Poly>
  Poly solve_impl(const array_const_view<double, 2, row_major> &rhs) const;

  int n_dims() const;

private:
  std::shared_ptr<Grid> grid;
  int_t i_cell;
  int order;

  LDLT ldlt;
  Eigen::MatrixXd A;
};

Eigen::MatrixXd assemble_weno_ao_matrix(const Grid &grid,
                                        const Stencil &stencil);

Eigen::MatrixXd assemble_weno_ao_matrix(
    const Grid &grid, const array_const_view<int_t, 1> &stencil, int order);

} // namespace zisa
#endif /* end of include guard */
