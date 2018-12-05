/* Generate TimeKeeper from CLI input parameters.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-03-07
 */
#ifndef TIME_KEEPER_INTERFACE_H_0KERGNWD
#define TIME_KEEPER_INTERFACE_H_0KERGNWD

#include <zisa/cli/input_parameters.hpp>
#include <zisa/ode/plotting_steps.hpp>
#include <zisa/ode/time_keeper.hpp>

namespace zisa {

struct PlottingStepsParameters {
  double fps;
  int_t steps_per_frame;

public:
  PlottingStepsParameters(const InputParameters &params);
};

struct TimeKeeperParameters {
  std::string wall_clock_time;
  double final_time;
  int_t total_steps;

public:
  TimeKeeperParameters(const InputParameters &params);
};

std::shared_ptr<PlottingSteps>
make_plotting_steps(const PlottingStepsParameters &plotting_params,
                    const TimeKeeperParameters &time_keeper_params);

std::shared_ptr<TimeKeeper>
make_time_keeper(const TimeKeeperParameters &params);

} // namespace zisa
#endif /* end of include guard: TIME_KEEPER_INTERFACE_H_0KERGNWD */
