/*
 *
 */

#ifndef WENO_AO_H_A1UM2
#define WENO_AO_H_A1UM2

#include <zisa/config.hpp>

#include <zisa/reconstruction/lsq_solver_family.hpp>
#include <zisa/reconstruction/stencil_family.hpp>
#include <zisa/reconstruction/weno_ao_params.hpp>

namespace zisa {

class WENO_AO {
private:
  static constexpr int MAX_DEGREE = 4;

public:
  WENO_AO(const std::shared_ptr<Grid> &grid,
          int_t i_cell,
          const WENO_AO_Params &params);

  Poly2D<MAX_DEGREE> reconstruct(const array<double, 1> &qbar) const;
  const std::vector<int_t> &local2global() const;

private:
  StencilFamily stencils;
  LSQSolverFamily lsq_solvers;

  std::vector<double> linear_weights;
  mutable array<double, 1> rhs;
};

// Eigen::MatrixXd assemble_weno_ao_matrix(const Grid &grid,
//                                         const array<LocalIndex, 1> &stencil,
//                                         int order,
//                                         double factor);

} // namespace zisa
#endif /* end of include guard */
