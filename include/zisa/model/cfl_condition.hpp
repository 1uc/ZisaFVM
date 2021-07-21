// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef CFL_CONDITION_DECL_H_OKTP9
#define CFL_CONDITION_DECL_H_OKTP9

#include <zisa/model/all_variables_fwd.hpp>

namespace zisa {

/// Interface for CFL-type time-step selection.
class CFLCondition {
public:
  virtual ~CFLCondition() = default;

  /// Compute largest stable time-step, due to CFL.
  virtual double operator()(const AllVariables &u) = 0;
};

} // namespace zisa

#endif
