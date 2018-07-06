/* Runtime selection of time integrator.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2015-11-18
 */
#include <zisa/ode/runge_kutta.hpp>
#include <zisa/ode/time_integration_factory.hpp>

namespace zisa {

std::shared_ptr<TimeIntegration>
make_time_integration(const std::string &desc,
                      const std::shared_ptr<RateOfChange> &rate_of_change,
                      const std::shared_ptr<BoundaryCondition> &bc,
                      const AllVariablesDimensions &dims) {
  if (desc.compare("ForwardEuler") == 0) {
    return std::make_shared<ForwardEuler>(rate_of_change, bc, dims);
  } else if (desc.compare("SSP2") == 0) {
    return std::make_shared<SSP2>(rate_of_change, bc, dims);
  } else if (desc.compare("SSP3") == 0) {
    return std::make_shared<SSP3>(rate_of_change, bc, dims);
  } else if (desc.compare("Wicker") == 0) {
    return std::make_shared<Wicker>(rate_of_change, bc, dims);
  } else if (desc.compare("RK4") == 0) {
    return std::make_shared<RK4>(rate_of_change, bc, dims);
  } else if (desc.compare("Fehlberg") == 0) {
    return std::make_shared<Fehlberg>(rate_of_change, bc, dims);
  }

  LOG_ERR("Unknown time integrator.");
}

} // namespace zisa
