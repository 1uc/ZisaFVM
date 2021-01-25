#include <zisa/io/scattered_data_source.hpp>

namespace zisa {

ScatteredDataSource::ScatteredDataSource(
    std::unique_ptr<AllVariablesScatterer> all_vars_scatterer,
    std::shared_ptr<Permutation> permutation,
    std::shared_ptr<DataSource> data_source,
    std::shared_ptr<HaloExchange> halo_exchange,
    const AllVariablesDimensions &all_var_dims)
    : scatterer(std::move(all_vars_scatterer)),
      permutation(std::move(permutation)),
      data_source(std::move(data_source)),
      halo_exchange(std::move(halo_exchange)) {

  LOG_ERR_IF(this->scatterer == nullptr, "Received a `nullptr`.");

  if (this->scatterer->is_this_rank_scattering()) {
    LOG_ERR_IF(this->permutation == nullptr, "Received a `nullptr`.");
    LOG_ERR_IF(this->data_source == nullptr, "Received a `nullptr`.");

    buffer = AllVariables(all_var_dims);
  }
}

void ScatteredDataSource::do_load(AllVariables &all_variables,
                                  SimulationClock &simulation_clock) {
  if (scatterer->is_this_rank_scattering()) {
    (*data_source)(buffer, simulation_clock);
    reverse_permutation(array_view(buffer.cvars), *permutation);
    reverse_permutation(array_view(buffer.avars), *permutation);

    scatterer->copy_local_patch(all_variables, buffer);
    scatterer->send(buffer);
  } else {
    scatterer->receive(all_variables);
  }

  (*halo_exchange)(all_variables);
  (*halo_exchange).wait();
}
}