#ifndef ZISA_DOMAIN_DECOMPOSITION_HPP_IICOX
#define ZISA_DOMAIN_DECOMPOSITION_HPP_IICOX

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/parallelization/mpi_halo_exchange.hpp>
#include <zisa/reconstruction/stencil_family.hpp>

namespace zisa {
struct PartitionedGrid {
  array<int_t, 1> partition;
  array<int_t, 1> boundaries;
  array<int_t, 1> permutation;
};

std::tuple<array<int_t, 2>, array<XYZ, 1>, array<StencilFamily, 1>, Halo>
extract_subgrid(const Grid &grid,
                const PartitionedGrid &partitioned_grid,
                const array<StencilFamily, 1> &stencils,
                int_t k_part);

}

#endif // ZISA_DOMAIN_DECOMPOSITION_HPP
