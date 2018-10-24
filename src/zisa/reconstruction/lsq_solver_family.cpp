#include <zisa/reconstruction/lsq_solver_family.hpp>

namespace zisa {

LSQSolverFamily::LSQSolverFamily(const std::shared_ptr<Grid> &grid,
                                 const StencilFamily &stencils) {

  assert(stencils.size() > 0);

  for (const auto &s : stencils) {
    solvers_.emplace_back(grid, s);
  }

  assert(solvers_.size() == stencils.size());
}

} // namespace zisa
