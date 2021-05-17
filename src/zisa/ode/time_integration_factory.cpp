#include <zisa/ode/runge_kutta.hpp>
#include <zisa/ode/time_integration_factory.hpp>

namespace zisa {

std::shared_ptr<TimeIntegration>
make_time_integration(const std::string &desc,
                      const std::shared_ptr<RateOfChange> &rate_of_change,
                      const std::shared_ptr<BoundaryCondition> &bc,
                      const AllVariablesDimensions &dims) {
  if (desc == "ForwardEuler") {
    return std::make_shared<ForwardEuler>(rate_of_change, bc, dims);
  } else if (desc == "SSP2") {
    return std::make_shared<SSP2>(rate_of_change, bc, dims);
  } else if (desc == "SSP3") {
    return std::make_shared<SSP3>(rate_of_change, bc, dims);
  } else if (desc == "Wicker") {
    return std::make_shared<Wicker>(rate_of_change, bc, dims);
  } else if (desc == "RK4") {
    return std::make_shared<RK4>(rate_of_change, bc, dims);
  } else if (desc == "Fehlberg") {
    return std::make_shared<Fehlberg>(rate_of_change, bc, dims);
  }

  LOG_ERR("Unknown time integrator.");
}

} // namespace zisa
