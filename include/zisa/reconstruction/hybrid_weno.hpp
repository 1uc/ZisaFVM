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

  array<double, 1> linear_weights;
  mutable array<double, 1> non_linear_weights;
  double epsilon;
  double exponent;

public:
  HybridWENO() = default;

  HybridWENO(const std::shared_ptr<Grid> &grid,
             int_t i_cell,
             const HybridWENOParams &params);

  HybridWENO(const std::shared_ptr<Grid> &grid,
             StencilFamily stencil_family,
             int_t i_cell,
             const HybridWENOParams &params);

  auto local2global() const -> decltype(stencils.local2global());
  int_t combined_stencil_size() const;

  /// Indistinguishable by calls to the public interface.
  bool operator==(const HybridWENO &other) const;

  /// Can be distinguished by calls to the public interface.
  bool operator!=(const HybridWENO &other) const;

  std::string str(int /* verbose */ = 0) const {
    return string_format("%s, eps = %e, s = %f",
                         format_as_list(linear_weights).c_str(),
                         epsilon,
                         exponent);
  }

protected:
  void compute_polys(array<double, 2, row_major> &rhs,
                     array<WENOPoly, 1> &polys,
                     const array<cvars_t, 1> &qbar) const;
  WENOPoly hybridize(array<WENOPoly, 1> &polys) const;
  WENOPoly eno_hybridize(array<WENOPoly, 1> &polys) const;
  WENOPoly tau_hybridize(array<WENOPoly, 1> &polys) const;
};

} // namespace zisa
#endif /* end of include guard */
