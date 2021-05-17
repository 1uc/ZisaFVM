#ifndef TIME_INTEGRATION_FACTORY_H_WA4FHC7U
#define TIME_INTEGRATION_FACTORY_H_WA4FHC7U

#include <zisa/config.hpp>

#include <zisa/boundary/boundary_condition.hpp>
#include <zisa/ode/rate_of_change.hpp>
#include <zisa/ode/runge_kutta.hpp>

namespace zisa {

std::shared_ptr<TimeIntegration>
make_time_integration(const std::string &desc,
                      const std::shared_ptr<RateOfChange> &rate_of_change,
                      const std::shared_ptr<BoundaryCondition> &bc,
                      const AllVariablesDimensions &dims);

} // namespace zisa

#endif /* end of include guard: TIME_INTEGRATION_FACTORY_H_WA4FHC7U */
