// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_INSTANTANEOUS_PHYSICS_NCJWO_HPP
#define ZISA_INSTANTANEOUS_PHYSICS_NCJWO_HPP

#include <zisa/model/all_variables_fwd.hpp>
#include <zisa/model/poisson_solver.hpp>
#include <zisa/ode/simulation_clock.hpp>

namespace zisa {
class InstantaneousPhysics {
public:
  virtual ~InstantaneousPhysics() = default;

  virtual void compute(const SimulationClock &clock,
                       AllVariables &all_variables)
      = 0;
};

class NoInstantaneousPhysics : public InstantaneousPhysics {
public:
  virtual void compute(const SimulationClock &, AllVariables &) override;
};

template <class EULER>
class SelfGravity : public InstantaneousPhysics {
private:
  using euler_t = EULER;

public:
  SelfGravity(std::shared_ptr<euler_t> euler,
              std::shared_ptr<PoissonSolver> poisson_solver)
      : euler(std::move(euler)), poisson_solver(std::move(poisson_solver)) {}

  void compute(const SimulationClock &, AllVariables &all_variables) override {
    poisson_solver->update(euler->gravity, all_variables);
  }

private:
  std::shared_ptr<euler_t> euler;
  std::shared_ptr<PoissonSolver> poisson_solver;
};

}

#endif // ZISA_INSTANEOUS_PHYSICS_HPP
