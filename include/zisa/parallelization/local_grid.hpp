#ifndef ZISA_LOCAL_GRID_HPP_UWDHT
#define ZISA_LOCAL_GRID_HPP_UWDHT

#include <zisa/grid/grid.hpp>
#include <zisa/parallelization/distributed_grid.hpp>
#include <zisa/parallelization/domain_decomposition.hpp>
#include <zisa/reconstruction/stencil_family.hpp>

namespace zisa {

std::tuple<std::shared_ptr<array<StencilFamily, 1>>,
           std::shared_ptr<DistributedGrid>,
           std::shared_ptr<Grid>>
load_local_grid(const std::string &subgrid_name,
                const StencilFamilyParams &stencil_params,
                const QRDegrees &qr_degrees,
                int mpi_rank);

}

#endif // ZISA_LOCAL_GRID_HPP
