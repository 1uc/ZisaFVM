#ifndef ZISA_INSTANTANEOUS_PHYSICS_NCJWO_HPP
#define ZISA_INSTANTANEOUS_PHYSICS_NCJWO_HPP

#include <zisa/model/all_variables_fwd.hpp>
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

class SelfGravity : public InstantaneousPhysics {};

}

#endif // ZISA_INSTANEOUS_PHYSICS_HPP
