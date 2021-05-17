#ifndef TIME_KEEPER_INTERFACE_H_0KERGNWD
#define TIME_KEEPER_INTERFACE_H_0KERGNWD

#include <zisa/cli/input_parameters.hpp>
#include <zisa/ode/plotting_steps.hpp>
#include <zisa/ode/time_keeper.hpp>

namespace zisa {

struct PlottingStepsParameters {
  int_t k0;
  double t0;
  double fps;
  int_t steps_per_frame;
  int_t n_snapshots;

public:
  explicit PlottingStepsParameters(const InputParameters &params);
};

struct TimeKeeperParameters {
  std::string wall_clock_time;
  double final_time;
  int_t total_steps;

public:
  explicit TimeKeeperParameters(const InputParameters &params);
};

std::shared_ptr<PlottingSteps>
make_plotting_steps(const PlottingStepsParameters &plotting_params,
                    const TimeKeeperParameters &time_keeper_params);

std::shared_ptr<TimeKeeper>
make_time_keeper(const TimeKeeperParameters &params);

} // namespace zisa
#endif /* end of include guard: TIME_KEEPER_INTERFACE_H_0KERGNWD */
