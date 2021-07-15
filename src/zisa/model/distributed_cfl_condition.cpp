// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/model/distributed_cfl_condition.hpp>

namespace zisa {

DistributedCFLCondition::DistributedCFLCondition(
    std::shared_ptr<CFLCondition> cfl_condition,
    std::shared_ptr<AllReduce> all_reduce)
    : cfl_condition(std::move(cfl_condition)),
      all_reduce(std::move(all_reduce)) {}

double DistributedCFLCondition::operator()(const AllVariables &u) {
  double dt_local = (*cfl_condition)(u);
  return (*all_reduce)(dt_local);
}

}