#ifndef LSQ_SOLVER_FAMILY_H_ORZ2M
#define LSQ_SOLVER_FAMILY_H_ORZ2M

#include <cassert>
#include <vector>

#include <zisa/config.hpp>
#include <zisa/reconstruction/lsq_solver.hpp>
#include <zisa/reconstruction/stencil_family.hpp>

namespace zisa {

class LSQSolverFamily {
public:
  LSQSolverFamily(const std::shared_ptr<Grid> &grid,
                  const StencilFamily &stencils);

  /// Returns the k-th LSQ solver.
  inline const LSQSolver &operator[](int_t k) const {
    assert(k < solvers_.size());
    return solvers_[k];
  }

  /// Returns the number of LSQ solvers.
  inline int_t size() const { return solvers_.size(); }

private:
  std::vector<LSQSolver> solvers_;
};

} // namespace zisa
#endif /* end of include guard */
