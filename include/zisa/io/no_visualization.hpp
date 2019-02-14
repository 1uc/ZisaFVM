#ifndef NO_VISUALIZATION_H_ELUV9
#define NO_VISUALIZATION_H_ELUV9

#include <zisa/io/visualization.hpp>

namespace zisa {

class NoVisualization : public Visualization {
protected:
  virtual void
  do_visualization(const AllVariables &all_variables,
                   const SimulationClock &simulation_clock) override;
};

} // namespace zisa
#endif
