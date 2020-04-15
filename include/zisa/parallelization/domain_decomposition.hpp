#ifndef ZISA_DOMAIN_DECOMPOSITION_HPP_IICOX
#define ZISA_DOMAIN_DECOMPOSITION_HPP_IICOX

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/parallelization/distributed_grid.hpp>
#include <zisa/parallelization/mpi_halo_exchange.hpp>
#include <zisa/reconstruction/stencil_family.hpp>

namespace zisa {
struct PartitionedGrid {
  PartitionedGrid(array<int_t, 1> partition,
                  array<int_t, 1> boundaries,
                  array<int_t, 1> permutation);

  array<int_t, 1> partition;
  array<int_t, 1> boundaries;
  array<int_t, 1> permutation;
};

std::map<int_t, int_t>
sparse_inverse_permutation(const array_const_view<int_t, 1> &sigma);

PartitionedGrid compute_partitioned_grid(
    const Grid &grid, const array<StencilFamily, 1> &stencils, int_t n_parts);

array<int_t, 2> renumbered_vertex_indices(const array<int_t, 2> &vertex_indices,
                                          const array<int_t, 1> &permutation);

std::tuple<array<int_t, 2>, array<XYZ, 1>, array<int_t, 1>>
extract_subgrid(const Grid &grid,
                const PartitionedGrid &partitioned_grid,
                const array<StencilFamily, 1> &stencils,
                int_t k_part);

std::tuple<array<int_t, 2>, array<XYZ, 1>, array<int_t, 1>>
extract_subgrid_v2(const Grid &grid,
                   const std::function<bool(int_t)> &is_inside);

template <class T>
array<T, 1> extract_subarray(const array<T, 1> &global_array,
                             const array<int_t, 1> &global_indices) {
  auto n_cells_local = global_indices.size();

  auto local_array = array<T, 1>(n_cells_local);
  for_each(index_range(n_cells_local),
           [&](int_t i) { local_array(i) = global_array(global_indices(i)); });

  return local_array;
}

array<StencilFamily, 1>
extract_stencils(const array<StencilFamily, 1> &global_stencils,
                 const array<int_t, 2> &global_neighbours,
                 const std::function<bool(int_t)> &is_inside,
                 const array<int_t, 1> &global_indices);

DistributedGrid
extract_distributed_subgrid(const DistributedGrid &dgrid,
                            const array<int_t, 1> &global_indices);

class StencilBasedIndicator {
public:
  StencilBasedIndicator(const Grid &grid,
                        const array<StencilFamily, 1> &stencils,
                        const std::function<bool(int_t)> &is_interior);

  bool operator()(int_t i) const;

protected:
  void append(const std::vector<int_t> &l2g);

private:
  std::vector<int_t> good_indices;
};

}
#endif // ZISA_DOMAIN_DECOMPOSITION_HPP
