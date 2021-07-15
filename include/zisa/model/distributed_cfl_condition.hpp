// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_DISTRIBUTED_CFL_CONDITION_HPP_CQPSI
#define ZISA_DISTRIBUTED_CFL_CONDITION_HPP_CQPSI

#include <zisa/config.hpp>
#include <zisa/model/cfl_condition.hpp>
#include <zisa/parallelization/all_reduce.hpp>

namespace zisa {

/// Compute on each part of the domain and then combine.
class DistributedCFLCondition : public CFLCondition {
public:
  DistributedCFLCondition(std::shared_ptr<CFLCondition> cfl_condition,
                          std::shared_ptr<AllReduce> all_reduce);

  virtual double operator()(const AllVariables &u) override;

private:
  std::shared_ptr<CFLCondition> cfl_condition;
  std::shared_ptr<AllReduce> all_reduce;
};

}

#endif // ZISA_DISTRIBUTED_CFL_CONDITION_HPP
