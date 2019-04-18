#include <zisa/model/instantaneous_physics.hpp>

namespace zisa {
void NoInstantaneousPhysics::compute(const zisa::SimulationClock &,
                                     zisa::AllVariables &) {
  return;
}
}