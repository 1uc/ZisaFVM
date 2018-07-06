/* Some general boundary conditions, e.g. outflow, reflective, etc.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2014-12-03
 */
#ifndef BOUNDARY_CONDITION_H_KLH67JU5
#define BOUNDARY_CONDITION_H_KLH67JU5

#include <zisa/config.hpp>
#include <zisa/model/all_variables_fwd.hpp>

namespace zisa {
/// Interface for applying ghost-cell boundary conditions.
class BoundaryCondition {
public:
  virtual ~BoundaryCondition() = default;

  /// Apply the boundary conditions to `u`.
  virtual void apply(AllVariables &u, double t) = 0;

  /// Self-documenting string.
  virtual std::string str(int indent = 0) const = 0;
};

} // namespace zisa

#endif /* end of include guard: BOUNDARY_CONDITION_H_KLH67JU5 */
