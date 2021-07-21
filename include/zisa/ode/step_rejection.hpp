// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_STEP_REJECTION_DIOWK_HPP
#define ZISA_STEP_REJECTION_DIOWK_HPP

#include <zisa/model/all_variables.hpp>

namespace zisa {

class StepRejection {
public:
  virtual ~StepRejection() = default;

  virtual bool is_good() const = 0;
  virtual void check(const AllVariables &u0, const AllVariables &u1) const = 0;
  virtual double pick_time_step(double dt) const = 0;
};

class RejectNothing : public StepRejection {
public:
  virtual bool is_good() const override;
  virtual void check(const AllVariables &, const AllVariables &) const override;

  virtual double pick_time_step(double dt) const override;
};

class RejectLargeDensityChange : public StepRejection {
public:
  explicit RejectLargeDensityChange(double drho_crit_rel, int max_exponent);

  virtual bool is_good() const override;

  virtual void check(const AllVariables &all_vars0,
                     const AllVariables &all_vars1) const override;

  virtual double pick_time_step(double dt) const override;

private:
  mutable bool is_good_ = true;

  double growth_factor = 2.0;
  mutable double current_factor = 1.0;

  double drho_crit_rel;
  double minimum_factor;
};

}

#endif // ZISA_STEP_REJECTION_HPP
