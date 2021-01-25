#include <zisa/parallelization/local_grid.hpp>

namespace zisa {

std::tuple<std::shared_ptr<array<StencilFamily, 1>>,
           std::shared_ptr<DistributedGrid>,
           std::shared_ptr<Grid>>
load_local_grid(const std::string &subgrid_name,
                const std::function<bool(const Grid &, int_t)> &boundary_mask,
                const StencilFamilyParams &stencil_params,
                const QRDegrees &qr_degrees,
                int mpi_rank) {
  auto super_subgrid = zisa::load_grid(subgrid_name, qr_degrees);
  auto super_sub_dgrid = zisa::load_distributed_grid(subgrid_name);

  auto is_interior
      = [mpi_rank, &partition = super_sub_dgrid.partition](int_t i) {
          return partition[i] == integer_cast<int_t>(mpi_rank);
        };

  mask_ghost_cells(*super_subgrid,
                   [&boundary_mask, &is_interior](const Grid &grid, int_t i) {
                     return boundary_mask(grid, i) || !is_interior(i);
                   });

  auto super_sub_stencils
      = compute_stencil_families(*super_subgrid, stencil_params);

  auto is_needed
      = StencilBasedIndicator(*super_subgrid, super_sub_stencils, is_interior);

  auto [local_vertex_indices, local_vertices, super_sub_indices]
      = extract_subgrid_v2(*super_subgrid, is_needed);

  auto stencils = std::make_shared<array<StencilFamily, 1>>(
      extract_stencils(super_sub_stencils,
                       super_subgrid->neighbours,
                       is_interior,
                       super_sub_indices));

  auto dgrid = std::make_shared<DistributedGrid>(
      extract_distributed_subgrid(super_sub_dgrid, super_sub_indices));

  auto element_type = super_subgrid->element_type();
  auto grid = std::make_shared<Grid>(element_type,
                                     std::move(local_vertices),
                                     std::move(local_vertex_indices),
                                     qr_degrees);

  return std::tuple{stencils, dgrid, grid};
}

}
