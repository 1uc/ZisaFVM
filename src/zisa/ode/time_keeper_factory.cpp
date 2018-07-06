/* Generate TimeKeeper from CLI input parameters.
 */
#include <zisa/ode/time_keeper.hpp>
#include <zisa/ode/time_keeper_factory.hpp>

namespace zisa {
std::shared_ptr<PlottingSteps>
make_plotting_steps(const TimeKeeperParameters &params, double t_zero) {

  auto &plot_params = params.plotting_steps_params;
  if (plot_params.mode == PlottingMode::fixed_time_steps) {
    return std::make_shared<PlotAtFixedTimeSteps>(plot_params.frames);
  } else if (plot_params.mode == PlottingMode::fixed_interval) {

    double t_end = params.use_fixed_duration ? params.final_time : 1e+300;
    return std::make_shared<PlotAtFixedInterval>(
        t_zero, 1.0 / plot_params.fps, t_end);
  } else if (plot_params.mode == PlottingMode::every_nth_step) {
    return std::make_shared<PlotEveryNthStep>(plot_params.steps_per_frame);
  } else {
    LOG_ERR("Unrecognised mode.")
  }
}

std::shared_ptr<TimeKeeper>
make_time_keeper(const TimeKeeperParameters &params) {
  if (!params.wall_clock_time.empty()) {
    if (params.use_fixed_duration) {
      return std::make_shared<FixedWallClock<FixedDuration>>(
          params.wall_clock_time, params.final_time);
    } else {
      return std::make_shared<FixedWallClock<FixedTimeSteps>>(
          params.wall_clock_time, params.total_steps);
    }
  } else {
    if (params.use_fixed_duration) {
      return std::make_shared<FixedDuration>(params.final_time);
    } else {
      return std::make_shared<FixedTimeSteps>(params.total_steps);
    }
  }
}

} // namespace zisa
