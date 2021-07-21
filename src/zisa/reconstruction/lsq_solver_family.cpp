// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

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

auto LSQSolverFamily::begin() -> decltype(solvers_.begin()) {
  return solvers_.begin();
}
auto LSQSolverFamily::begin() const -> decltype(solvers_.begin()) {
  return solvers_.begin();
}

auto LSQSolverFamily::end() -> decltype(solvers_.end()) {
  return solvers_.end();
}
auto LSQSolverFamily::end() const -> decltype(solvers_.end()) {
  return solvers_.end();
}

bool operator==(const LSQSolverFamily &lhs, const LSQSolverFamily &rhs) {
  return std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

bool operator!=(const LSQSolverFamily &lhs, const LSQSolverFamily &rhs) {
  return !(lhs == rhs);
}

} // namespace zisa
