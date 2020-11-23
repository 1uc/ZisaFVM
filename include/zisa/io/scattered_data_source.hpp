#ifndef ZISA_SCATTERED_DATA_SOURCE_HPP_VIOEW
#define ZISA_SCATTERED_DATA_SOURCE_HPP_VIOEW

#include <zisa/config.hpp>

#include <thread>
#include <zisa/io/data_source.hpp>
#include <zisa/math/permutation.hpp>
#include <zisa/parallelization/all_variables_scatterer.hpp>
#include <zisa/parallelization/halo_exchange.hpp>

namespace zisa {

/// Sources data only on a subset of ranks.
/** The idea is to scatter the data obtained on specific ranks from those ranks.
 */
class ScatteredDataSource : public DataSource {
public:
  ScatteredDataSource(std::unique_ptr<AllVariablesScatterer> all_vars_scatterer,
                      std::shared_ptr<Permutation> permutation,
                      std::shared_ptr<DataSource> data_source,
                      std::shared_ptr<HaloExchange> halo_exchange,
                      const AllVariablesDimensions &all_var_dims);

protected:
  void do_load(AllVariables &all_variables,
               SimulationClock &simulation_clock) override;

private:
  AllVariables buffer;
  std::unique_ptr<AllVariablesScatterer> scatterer;
  std::shared_ptr<Permutation> permutation;
  std::shared_ptr<DataSource> data_source;
  std::shared_ptr<HaloExchange> halo_exchange;
};

}

#endif // ZISA_SCATTERED_DATA_SOURCE_HPP
