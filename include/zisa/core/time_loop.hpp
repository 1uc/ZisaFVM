/* Finite Volume Methods using CUDA.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2014-12-11
 */
#ifndef TIME_LOOP_H_3IJELTQK
#define TIME_LOOP_H_3IJELTQK

#include <zisa/model/sanity_check.hpp>
#include <zisa/datetime.hpp>
#include <zisa/io/progress_bar.hpp>
#include <zisa/io/visualization.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/ode/simulation_clock.hpp>
#include <zisa/ode/time_integration.hpp>
#include <zisa/model/cfl_condition.hpp>

namespace zisa {
/// The time-stepping loop.
/** Advance a solution forwards in time by repeatedly calling a
 *  one step solver.
 */
class TimeLoop {
public:
  TimeLoop(const std::shared_ptr<TimeIntegration> &time_integration,
           const std::shared_ptr<SimulationClock> &simulation_clock,
           const std::shared_ptr<CFLCondition> &cfl_condition,
           const std::shared_ptr<SanityCheck> &sanity_check,
           const std::shared_ptr<Visualization> &visualization);

  virtual ~TimeLoop() = default;

  /// Advance `u0` forwards in time.
  /** @param u0  Initial conditions in host memory.
   */
  void operator()(std::shared_ptr<AllVariables> u0);

  /// Things to do before entering the time loop.
  virtual void pre_loop(AllVariables &) {
    // empty default
  }

  /// Things to do before the update step.
  virtual void pre_update(AllVariables &) {
    // empty default
  }

  /// Things to do just after 'u0' has been updated.
  virtual void post_update(AllVariables &u0);

  /// Things to do after the time loop.
  virtual void post_loop(AllVariables &) {
    // empty default
  }

  /// Self-documenting string.
  virtual std::string str() const;

protected:
  void write_output(const AllVariables &u0);

  virtual void print_welcome_message(void) const;
  virtual void print_progress_message(void);
  virtual void print_goodbye_message(void) const;

  virtual void sanity_check(const AllVariables &all_variables) const;
  virtual void pick_time_step(const AllVariables &all_variables);

  void start_timer();
  void stop_timer();

protected:
  std::shared_ptr<TimeIntegration> time_integration;
  std::shared_ptr<SimulationClock> simulation_clock;
  std::shared_ptr<Visualization> visualization;
  std::shared_ptr<CFLCondition> cfl_condition;
  std::shared_ptr<SanityCheck> is_sane;

  ProgressBar progress_bar;

  time_stamp_t start_time;
  time_stamp_t end_time;
};

} // namespace zisa
#endif /* end of include guard */