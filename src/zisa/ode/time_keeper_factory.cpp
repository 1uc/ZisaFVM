/* Generate TimeKeeper from CLI input parameters.
 */
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/ode/time_keeper.hpp>
#include <zisa/ode/time_keeper_factory.hpp>

namespace zisa {
PlottingStepsParameters::PlottingStepsParameters(const InputParameters &params)
    : k0(0),
      t0(0.0),
      fps(-1.0),
      steps_per_frame(int_t(-1)),
      n_snapshots(int_t(-1)) {

  LOG_ERR_IF(!has_key(params, "io"), "Missing config section 'io'.");

  const auto &params_plotting = params["io"];

  //  if (has_key(params, "restart")) {
  //    auto filename = std::string(params["restart"]["file"]);
  //    auto reader = HDF5SerialReader(filename);
  //
  //    t0 = reader.read_scalar<double>("time");
  //    k0 = reader.read_scalar<int_t>("n_steps");
  //  }

  if (has_key(params_plotting, "fps")) {
    fps = params_plotting["fps"];
  }

  if (has_key(params_plotting, "steps_per_frame")) {
    steps_per_frame = params_plotting["steps_per_frame"];
  }

  if (has_key(params_plotting, "n_snapshots")) {
    n_snapshots = params_plotting["n_snapshots"];
  }
}

std::shared_ptr<PlottingSteps>
make_plotting_steps(const PlottingStepsParameters &plotting_params,
                    const TimeKeeperParameters &time_keeper_params) {

  double t_end = time_keeper_params.final_time;

  double fps = plotting_params.fps;
  int_t steps_per_frame = plotting_params.steps_per_frame;
  int_t n_snapshots = plotting_params.n_snapshots;

  if (fps > 0.0 && std::isfinite(t_end)) {
    double t0 = plotting_params.t0;
    double dt = 1.0 / fps;

    return std::make_shared<PlotAtFixedInterval>(t0, dt, t_end);
  }

  if (n_snapshots != int_t(-1) && std::isfinite(t_end)) {
    double t0 = plotting_params.t0;
    double dt = t_end / double(n_snapshots);

    return std::make_shared<PlotAtFixedInterval>(t0, dt, t_end);
  }

  if (steps_per_frame != int_t(-1) /* which might be > 0 */) {
    int_t k0 = plotting_params.k0;
    return std::make_shared<PlotEveryNthStep>(k0, steps_per_frame);
  }

  LOG_ERR("Could not decide on a plotting mode.");
}

TimeKeeperParameters::TimeKeeperParameters(const InputParameters &params)
    : final_time(std::numeric_limits<double>::max()),
      total_steps(std::numeric_limits<int_t>::max()) {

  LOG_ERR_IF(!has_key(params, "time"), "Failed to find section 'time'.");

  const auto &params_time = params["time"];

  final_time = std::numeric_limits<double>::infinity();
  for (auto &&key : std::vector<std::string>{"t_end", "final_time"}) {
    if (has_key(params_time, key)) {
      final_time = params_time[key];
    }
  }

  total_steps = std::numeric_limits<int_t>::max();
  if (has_key(params_time, "n_steps")) {
    total_steps = params_time["n_steps"];
  }

  if (has_key(params_time, "wall_clock")) {
    wall_clock_time = params_time["wall_clock"];
  }
}

std::shared_ptr<TimeKeeper>
make_time_keeper(const TimeKeeperParameters &params) {
  auto final_time = params.final_time;
  auto total_steps = params.total_steps;
  auto wall_clock = params.wall_clock_time;

  if (!wall_clock.empty()) {
    using time_keeper_t = FixedWallClock<FixedDurationAndTimeSteps>;
    return std::make_shared<time_keeper_t>(wall_clock, final_time, total_steps);
  } else {
    using time_keeper_t = FixedDurationAndTimeSteps;
    return std::make_shared<time_keeper_t>(final_time, total_steps);
  }
}

} // namespace zisa
