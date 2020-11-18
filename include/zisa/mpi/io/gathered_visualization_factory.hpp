#ifndef ZISA_GATHERED_VISUALIZATION_FACTORY_HPP_EUYBY
#define ZISA_GATHERED_VISUALIZATION_FACTORY_HPP_EUYBY

#include <zisa/io/gathered_visualization.hpp>
#include <zisa/mpi/io/gathered_vis_info.hpp>
#include <zisa/mpi/parallelization/mpi_single_node_array_gatherer.hpp>

namespace zisa {

std::shared_ptr<GatheredVisualization> make_gathered_visualization(
    std::shared_ptr<GatheredVisInfo> vis_info,
    std::shared_ptr<MPISingleNodeArrayGathererFactory> gatherer_factory,
    std::shared_ptr<Visualization> visualization,
    AllVariablesDimensions all_var_dims);

}

#endif // ZISA_GATHERED_VISUALIZATION_FACTORY_HPP
