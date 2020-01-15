/* Finite Volume Methods using CUDA.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2014-12-11
 */
#ifndef TIME_LOOP_H_3IJELTQK
#define TIME_LOOP_H_3IJELTQK

#include <zisa/datetime.hpp>
#include <zisa/io/progress_bar.hpp>
#include <zisa/io/visualization.hpp>
#include <zisa/model/cfl_condition.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/instantaneous_physics.hpp>
#include <zisa/model/sanity_check.hpp>
#include <zisa/ode/simulation_clock.hpp>
#include <zisa/ode/step_rejection.hpp>
#include <zisa/ode/time_integration.hpp>

namespace zisa {
/// The time-stepping loop.
/** Advance a solution forwards in time by repeatedly calling a
 *  one step solver.
 */
class TimeLoop {
public:
  TimeLoop(const std::shared_ptr<TimeIntegration> &time_integration,
           const std::shared_ptr<InstantaneousPhysics> &instantaneous_physics,
           const std::shared_ptr<StepRejection> &step_rejection,
           const std::shared_ptr<SimulationClock> &simulation_clock,
           const std::shared_ptr<CFLCondition> &cfl_condition,
           const std::shared_ptr<SanityCheck> &sanity_check,
           const std::shared_ptr<Visualization> &visualization,
           const std::shared_ptr<ProgressBar> &progress_bar);

  virtual ~TimeLoop() = default;

  /// Advance `u0` forwards in time.
  /** @param u0  Initial conditions in host memory.
   */
  std::shared_ptr<AllVariables> operator()(std::shared_ptr<AllVariables> u0);

  /// Self-documenting string.
  virtual std::string str() const;

protected:
  virtual void reject_step(std::shared_ptr<AllVariables> &u0,
                           std::shared_ptr<AllVariables> &u1);

  void write_output(const AllVariables &u0);
  void post_update(AllVariables &u0);

  virtual void print_welcome_message() const;
  virtual void print_progress_message();
  virtual void print_goodbye_message() const;
  virtual void print(const std::string &str) const;

  virtual void sanity_check(const AllVariables &all_variables) const;
  virtual void pick_time_step(const AllVariables &all_variables);

  void start_timer();
  void stop_timer();

protected:
  std::shared_ptr<TimeIntegration> time_integration;
  std::shared_ptr<InstantaneousPhysics> instantaneous_physics;
  std::shared_ptr<StepRejection> step_rejection;
  std::shared_ptr<SimulationClock> simulation_clock;
  std::shared_ptr<Visualization> visualization;
  std::shared_ptr<CFLCondition> cfl_condition;
  std::shared_ptr<SanityCheck> is_sane;
  std::shared_ptr<ProgressBar> progress_bar;

  time_stamp_t start_time;
  time_stamp_t end_time;
};

} // namespace zisa
#endif /* end of include guard */
