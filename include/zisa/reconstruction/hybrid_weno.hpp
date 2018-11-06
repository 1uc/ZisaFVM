/* Base class for hybridized WENO.
 */

#ifndef HYBRID_WENO_H_49M5Q
#define HYBRID_WENO_H_49M5Q

#include <zisa/config.hpp>

#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/hybrid_weno_params.hpp>
#include <zisa/reconstruction/lsq_solver_family.hpp>
#include <zisa/reconstruction/stencil_family.hpp>

namespace zisa {

class HybridWENO {
private:
  static constexpr int MAX_DEGREE = 4;

protected:
  StencilFamily stencils;
  LSQSolverFamily lsq_solvers;

  mutable array<Poly2D<MAX_DEGREE>, 1> polys;
  mutable array<double, 1> rhs;

  array<double, 1> linear_weights;
  double epsilon;
  double exponent;

public:
  HybridWENO(const std::shared_ptr<Grid> &grid,
             int_t i_cell,
             const HybridWENO_Params &params);

  Poly2D<MAX_DEGREE> reconstruct(const array<double, 1> &qbar) const;
  auto local2global() const -> decltype(stencils.local2global());

  /// Indistinguishable by calls to the public interface.
  bool operator==(const HybridWENO &other) const;

  /// Can be distinguished by calls to the public interface.
  bool operator!=(const HybridWENO &other) const;

protected:
  void compute_polys(const array<double, 1> &qbar) const;
  auto hybridize() const -> Poly2D<MAX_DEGREE>;
};

} // namespace zisa
#endif /* end of include guard */
