#ifndef LSQ_SOLVER_H_7JQG4
#define LSQ_SOLVER_H_7JQG4

#include <Eigen/Dense>

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/stencil.hpp>

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
  using QR = Eigen::FullPivHouseholderQR<Eigen::MatrixXd>;

private:
  /// Highest degree polynomial possible.
  static constexpr int MAX_DEGREE = 4;

public:
  LSQSolver(const std::shared_ptr<Grid> &grid, const Stencil &stencil);

  /// Solve the LSQ problem with right-hand side `rhs`.
  Poly2D<MAX_DEGREE> solve(const array<double, 1> &rhs) const;

private:
  std::shared_ptr<Grid> grid;
  int_t i_cell;
  int order;

  QR qr;
};

} // namespace zisa
#endif /* end of include guard */
