/* Generate TimeKeeper from CLI input parameters.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-03-07
 */
#ifndef TIME_KEEPER_INTERFACE_H_0KERGNWD
#define TIME_KEEPER_INTERFACE_H_0KERGNWD
#include <zisa/ode/plotting_steps.hpp>
#include <zisa/ode/time_keeper.hpp>

namespace zisa {

struct PlottingStepsParameters {
  PlottingMode mode;

  std::vector<int_t> frames;
  int_t steps_per_frame;
  double fps;
  double starting_time = 0.0;
};

struct TimeKeeperParameters {
  bool use_fixed_duration;

  std::string wall_clock_time;
  double final_time;
  int_t total_steps;

  PlottingStepsParameters plotting_steps_params;
};

std::shared_ptr<PlottingSteps>
make_plotting_steps(const TimeKeeperParameters &params, double t_zero);

std::shared_ptr<TimeKeeper>
make_time_keeper(const TimeKeeperParameters &params);

} // namespace zisa
#endif /* end of include guard: TIME_KEEPER_INTERFACE_H_0KERGNWD */
