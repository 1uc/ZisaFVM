#ifndef ZISA_MPI_DUMP_SNAPSHOT_HPP_CUMJS
#define ZISA_MPI_DUMP_SNAPSHOT_HPP_CUMJS

#include <thread>
#include <zisa/config.hpp>
#include <zisa/io/visualization.hpp>
#include <zisa/parallelization/all_variables_gatherer.hpp>

namespace zisa {

/// Performs visualization only on a subset of ranks.
/** The idea is to gather the data on specific ranks, and then
 *  run the desired visualization on a separate thread on those ranks.
 *
 *  The logic of gathering the data is abstracted in `ArrayGatherer`.
 *  The logic for visualization is abstracted in `Visualization`.
 */
class GatheredVisualization : public Visualization {
public:
  GatheredVisualization(std::shared_ptr<AllVariablesGatherer> all_vars_gatherer,
                        std::shared_ptr<Visualization> visualization,
                        const AllVariablesDimensions &all_var_dims);

  ~GatheredVisualization() override;

protected:
  void do_visualization(const AllVariables &all_variables,
                        const SimulationClock &simulation_clock) override;

private:
  AllVariables buffer;
  std::shared_ptr<AllVariablesGatherer> gatherer;
  std::shared_ptr<Visualization> visualization;
  std::unique_ptr<std::thread> job = nullptr;
};

}
#endif // ZISA_MPI_DUMP_SNAPSHOT_HPP
