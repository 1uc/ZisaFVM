#ifndef ZISA_SCATTERED_DATA_SOURCE_FACTORY_HPP_BQPUY
#define ZISA_SCATTERED_DATA_SOURCE_FACTORY_HPP_BQPUY

#include <zisa/config.hpp>

#include <zisa/io/scattered_data_source.hpp>
#include <zisa/mpi/io/gathered_vis_info.hpp>
#include <zisa/mpi/parallelization/mpi_single_node_array_scatterer.hpp>

namespace zisa {

std::shared_ptr<ScatteredDataSource> make_scattered_data_source(
    std::shared_ptr<GatheredVisInfo> vis_info,
    std::shared_ptr<MPISingleNodeArrayScattererFactory> scatterer_factory,
    std::shared_ptr<DataSource> data_source,
    AllVariablesDimensions all_var_dims);

}


#endif // ZISA_SCATTERED_DATA_SOURCE_FACTORY_HPP
