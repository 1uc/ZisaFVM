/* Base class for hybridized WENO.
 */

#ifndef HYBRID_WENO_H_49M5Q
#define HYBRID_WENO_H_49M5Q

#include <zisa/config.hpp>

#include <zisa/math/poly2d.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/reconstruction/hybrid_weno_params.hpp>
#include <zisa/reconstruction/lsq_solver_family.hpp>
#include <zisa/reconstruction/stencil_family.hpp>
#include <zisa/reconstruction/weno_poly.hpp>

namespace zisa {
class HybridWENO {
protected:
  using cvars_t = euler_var_t;

protected:
  StencilFamily stencils;
  LSQSolverFamily lsq_solvers;

  mutable array<WENOPoly, 1> polys;
  mutable array<double, 2, row_major> rhs;

  array<double, 1> linear_weights;
  mutable array<double, 1> non_linear_weights;
  double epsilon;
  double exponent;

  int_t i_cell;

public:
  HybridWENO() = default;

  HybridWENO(const std::shared_ptr<Grid> &grid,
             int_t i_cell,
             const HybridWENOParams &params);

  auto local2global() const -> decltype(stencils.local2global());
  int_t combined_stencil_size() const;

  /// Indistinguishable by calls to the public interface.
  bool operator==(const HybridWENO &other) const;

  /// Can be distinguished by calls to the public interface.
  bool operator!=(const HybridWENO &other) const;

protected:
  void compute_polys(const array<cvars_t, 1> &qbar) const;
  WENOPoly hybridize() const;
  WENOPoly eno_hybridize() const;
  WENOPoly tau_hybridize() const;
};

} // namespace zisa
#endif /* end of include guard */
