/* Compiled versions of FVM.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-08-25
 */

#include <zisa/config.hpp>

#include <zisa/core/time_loop.hpp>

#include <zisa/math/cartesian.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/sanity_check.hpp>

namespace zisa {
TimeLoop::TimeLoop(const std::shared_ptr<TimeIntegration> &time_integration,
                   const std::shared_ptr<SimulationClock> &simulation_clock,
                   const std::shared_ptr<CFLCondition> &cfl_condition,
                   const std::shared_ptr<SanityCheck> &sanity_check,
                   const std::shared_ptr<Visualization> &visualization)
    : time_integration(time_integration),
      simulation_clock(simulation_clock),
      visualization(visualization),
      cfl_condition(cfl_condition),
      is_sane(sanity_check),
      progress_bar(1) {}

void TimeLoop::operator()(std::shared_ptr<AllVariables> u0) {
  start_timer();
  print_welcome_message();
  progress_bar.reset();

  write_output(*u0);
  sanity_check(*u0);

  pick_time_step(*u0);

  // pre loop routines, set up additional state, etc.
  pre_loop(*u0);

  while (!simulation_clock->is_finished()) {
    // pre update routines (mainly for increased flexibility).
    pre_update(*u0);

    double t = simulation_clock->current_time();
    double dt = simulation_clock->current_time_step();

    // v-- possibly different buffer, contains u1.
    u0 = time_integration->compute_step(u0, t, dt);
    simulation_clock->advance();
    pick_time_step(*u0);

    // post update routines, plotting, etc.
    post_update(*u0);

    print_progress_message();
  }

  // post loop routines, for things like start-end comparisons.
  post_loop(*u0);
  stop_timer();
  print_goodbye_message();
}

void TimeLoop::pick_time_step(const AllVariables &all_variables) {
  double dt_cfl = (*cfl_condition)(all_variables);
  simulation_clock->set_time_step(dt_cfl);
}

void TimeLoop::post_update(AllVariables &u0) {
  write_output(u0);
  sanity_check(u0);
}

void TimeLoop::print_welcome_message(void) const {
  std::string date = date_format(time(0));

  if (simulation_clock->current_time() == 0.0) {
    std::cout << "-------- New run --------- \n";
  } else {
    std::cout << "-------- Continued run --------- \n";
  }
  std::cout << "     Date: " << date << "\n";
  std::cout << "-------- Solver --------\n";
  std::cout << str() << "\n";
  std::cout << "------------------------\n";
}

void TimeLoop::print_progress_message(void) {
  progress_bar.write_progress(std::cout,
                              simulation_clock->compact_progess_string());
}

void TimeLoop::print_goodbye_message(void) const {
  std::cout << "\n ---- ----- ----\n";
  std::cout << "total steps: " << simulation_clock->current_step() << "\n";
  std::cout << " final time: " << simulation_clock->current_time() << "\n";

  std::string date = date_format(time(0));

  auto elapsed = elapsed_seconds(end_time, start_time);
  std::string duration = duration_format(elapsed);

  std::cout << "-------------------------- \n";
  std::cout << "     Date: " << date << "\n";
  std::cout << " Duration: " << duration << "\n";
  if (simulation_clock->is_interrupted()) {
    std::cout << "-------- to be continued --------- \n";
  } else {
    std::cout << "-------- End of run --------- \n";
  }
}

void TimeLoop::start_timer() { start_time = current_time_stamp(); }

void TimeLoop::stop_timer() { end_time = current_time_stamp(); }

void TimeLoop::write_output(const AllVariables &all_variables) {
  if (simulation_clock->is_plotting_step()) {
    (*visualization)(all_variables, *simulation_clock);
  }
}

void TimeLoop::sanity_check(const AllVariables &all_variables) const {
  if (!(*is_sane)(all_variables)) {
    (*visualization)(all_variables, *simulation_clock);
    LOG_ERR("Values aren't plausible.\n");
  }
}

std::string TimeLoop::str() const {
  return "FVM (`TimeLoop`)\n" + indent_block(1, time_integration->str());
}

} // namespace zisa
