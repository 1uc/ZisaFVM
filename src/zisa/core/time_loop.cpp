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

#if ZISA_HAS_MPI == 1
#include <zisa/mpi/mpi.hpp>
#endif

namespace zisa {

TimeLoop::TimeLoop(
    const std::shared_ptr<TimeIntegration> &time_integration,
    const std::shared_ptr<InstantaneousPhysics> &instantaneous_physics,
    const std::shared_ptr<StepRejection> &step_rejection,
    const std::shared_ptr<SimulationClock> &simulation_clock,
    const std::shared_ptr<CFLCondition> &cfl_condition,
    const std::shared_ptr<SanityCheck> &sanity_check,
    const std::shared_ptr<Visualization> &visualization,
    const std::shared_ptr<ProgressBar> &progress_bar)
    : time_integration(time_integration),
      instantaneous_physics(instantaneous_physics),
      step_rejection(step_rejection),
      simulation_clock(simulation_clock),
      visualization(visualization),
      cfl_condition(cfl_condition),
      is_sane(sanity_check),
      progress_bar(progress_bar) {}

std::shared_ptr<AllVariables>
TimeLoop::operator()(std::shared_ptr<AllVariables> u0) {
  start_timer();
  print_welcome_message();
  progress_bar->reset();

  write_output(*u0);
  sanity_check(*u0);

  pick_time_step(*u0);

  while (!simulation_clock->is_finished()) {
    double t = simulation_clock->current_time();
    double dt = simulation_clock->current_time_step();

    // v-- possibly different buffer, contains u1.
    auto u1 = time_integration->compute_step(u0, t, dt);

    step_rejection->check(*u0, *u1);
    if (!step_rejection->is_good()) {
      reject_step(u0, u1);
      continue;
    }

    instantaneous_physics->compute(*simulation_clock, *u1);

    simulation_clock->advance();
    pick_time_step(*u1);

    post_update(*u1);
    print_progress_message();

    u0 = u1;
  }

  stop_timer();
  print_goodbye_message();

  return u0;
}

void TimeLoop::pick_time_step(const AllVariables &all_variables) {
  double dt_cfl = (*cfl_condition)(all_variables);
  double dt_safe = step_rejection->pick_time_step(dt_cfl);
  simulation_clock->set_time_step(dt_safe);
}

void TimeLoop::post_update(AllVariables &u0) {
  write_output(u0);
  sanity_check(u0);
}

void TimeLoop::print_welcome_message() const {
  std::stringstream ss;

  std::string date = date_format(time(nullptr));
  if (simulation_clock->current_time() == 0.0) {
    ss << "-------- New run --------- \n";
  } else {
    ss << "-------- Continued run --------- \n";
  }
  ss << "     Date: " << date << "\n";
  ss << "-------- Solver --------\n";
  ss << str() << "\n";
  ss << "------------------------\n";

  print(ss.str());
}

void TimeLoop::print_progress_message() {
  progress_bar->write_progress(std::cout,
                               simulation_clock->compact_progess_string());
}

void TimeLoop::print_goodbye_message() const {
  std::stringstream ss;

  ss << "\n ---- ----- ----\n";
  ss << "total steps: " << simulation_clock->current_step() << "\n";
  ss << " final time: " << simulation_clock->current_time() << "\n";

  std::string date = date_format(time(nullptr));

  auto elapsed = elapsed_seconds(end_time, start_time);
  std::string duration = duration_format(elapsed);

  ss << "-------------------------- \n";
  ss << "     Date: " << date << "\n";
  ss << " Duration: " << duration << "\n";
  ss << "  Elapsed: " << elapsed << " s\n";
  if (simulation_clock->is_interrupted()) {
    ss << "-------- to be continued --------- \n";
  } else {
    ss << "-------- End of run --------- \n";
  }

  auto of = std::ofstream("run_time.txt");
  of << elapsed << "\n";

  print(ss.str());
}

void TimeLoop::print(const std::string &str) const {
#if ZISA_HAS_MPI == 1
  static int mpi_rank = zisa::mpi::rank(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << str;
  }
#else
  std::cout << str;
#endif
}

void TimeLoop::start_timer() {
  MPI_Barrier(MPI_COMM_WORLD); // FIXME
  start_time = current_time_stamp();
}

void TimeLoop::stop_timer() {
  MPI_Barrier(MPI_COMM_WORLD); // FIXME
  end_time = current_time_stamp();
}

void TimeLoop::write_output(const AllVariables &all_variables) {
  if (simulation_clock->is_plotting_step()) {
    (*visualization)(all_variables, *simulation_clock);
  }
}

void TimeLoop::sanity_check(const AllVariables &all_variables) const {
  if (!(*is_sane)(all_variables)) {
    //    (*visualization)(all_variables, *simulation_clock);
    //    visualization->wait();
    LOG_ERR(string_format("[%d] Values aren't plausible.\n",
                          simulation_clock->current_step()));
  }
}

std::string TimeLoop::str() const {
  return "FVM (`TimeLoop`)\n" + indent_block(1, time_integration->str());
}

void TimeLoop::reject_step(std::shared_ptr<AllVariables> &u0,
                           std::shared_ptr<AllVariables> &u1) {

  double dt = simulation_clock->current_time_step();
  auto dt_new = step_rejection->pick_time_step(dt);
  simulation_clock->set_time_step(dt_new);

  // Remember, `time_integration` keeps the buffer it
  // was passed for future use. Hence we must either copy data and copy the
  // pointers, or swap the buffers & copy the pointers.
  //
  // Note: That `u1` will not be used after this call.
  *u1 = *u0; // TODO implement swap(*u0, *u1)
  u0 = u1;
}

} // namespace zisa
