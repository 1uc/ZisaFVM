/* Interface for visualizing the flow.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-09-14
 */
#ifndef VISUALIZATION_H_748ODJNR
#define VISUALIZATION_H_748ODJNR

#include <zisa/config.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/ode/simulation_clock.hpp>

namespace zisa {
class Visualization {
public:
  virtual ~Visualization() = default;

  /// Perform visualization.
  void operator()(const AllVariables &all_variables,
                  const SimulationClock &simulation_clock);

  /// Visualize the steady-state.
  void steady_state(const AllVariables &all_variables);

  /// Wait until the visualization has completed.
  void wait();

protected:
  virtual void do_visualization(const AllVariables &all_variables,
                                const SimulationClock &simulation_clock)
      = 0;

  virtual void do_steady_state(const AllVariables &);

  virtual void do_wait();
};

} // namespace zisa
#endif /* end of include guard: VISUALIZATION_H_748ODJNR */
