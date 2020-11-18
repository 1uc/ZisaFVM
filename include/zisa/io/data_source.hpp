//
// Created by lucg on 11/17/20.
//

#ifndef ZISA_DATA_SOURCE_HPP
#define ZISA_DATA_SOURCE_HPP

#include <zisa/config.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/ode/simulation_clock.hpp>

namespace zisa {

class DataSource {
public:
  virtual ~DataSource() = default;

  /// Get the data from somewhere.
  /** This is also expected to set the simulation clock
   *  to the appropriate time and time-step.
   */
  void operator()(AllVariables &all_variables,
                  SimulationClock &simulation_clock);

protected:
  virtual void do_load(AllVariables &all_variables,
                       SimulationClock &simulation_clock)
      = 0;
};
}

#endif // ZISA_DATA_SOURCE_HPP
