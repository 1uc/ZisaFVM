#include <zisa/parallelization/domain_decomposition.hpp>

#include <map>
#include <metis.h>
#include <zisa/grid/neighbour_range.hpp>

namespace zisa {
using metis_idx_t = ::idx_t;

PartitionedGrid::PartitionedGrid(array<int_t, 1> partition,
                                 array<int_t, 1> boundaries,
                                 array<int_t, 1> permutation)
    : partition(std::move(partition)),
      boundaries(std::move(boundaries)),
      permutation(std::move(permutation)) {}

array<int_t, 1> compute_partition_full_stencil(
    const Grid &grid,
    const std::vector<std::vector<int_t>> & /*stencils*/,
    int n_parts) {

  auto n_cells = grid.n_cells;
  auto nvtxs = metis_idx_t(n_cells);
  auto ncon = metis_idx_t(1); // # constraints

  std::vector<std::vector<int_t>> graph(n_cells);
  for (int_t i = 0; i < n_cells; ++i) {
    // FIXME: currently nearest neighbours
    //    for (int_t k = 0; k < stencils[i].size(); ++k) {
    for (int_t k = 0; k < grid.max_neighbours; ++k) {
      if (!grid.is_valid(i, k)) {
        continue;
      }

      int_t j = grid.neighbours(i, k);
      // int_t j = stencils[i][k];

      if (std::find(graph[i].cbegin(), graph[i].cend(), j) == graph[i].cend()) {
        graph[i].push_back(j);
        graph[j].push_back(i);
      }
    }
  }

  int_t n_edges = 0;
  for (int_t i = 0; i < n_cells; ++i) {
    n_edges += graph[i].size();
    std::sort(graph[i].begin(), graph[i].end());
  }

  auto xadj = array<metis_idx_t, 1>(n_cells + 1);
  auto adjncy = array<metis_idx_t, 1>(n_edges);

  int_t k_edge = 0;
  for (int_t i = 0; i < n_cells; ++i) {
    xadj[i] = integer_cast<metis_idx_t>(k_edge);
    for (int_t j : graph[i]) {
      adjncy[k_edge] = integer_cast<metis_idx_t>(j);
      k_edge += 1;
    }
  }
  xadj[n_cells] = integer_cast<metis_idx_t>(n_edges);

  auto nparts = metis_idx_t(n_parts);
  auto objval = metis_idx_t(-1);

  array<metis_idx_t, 1> part(n_cells);

  // clang-format off
  METIS_PartGraphKway(
      &nvtxs, &ncon, xadj.raw(), adjncy.raw(), nullptr, nullptr, nullptr,
      &nparts, nullptr, nullptr, nullptr, &objval, part.raw());
  // clang-format on

  array<int_t, 1> partition(n_cells);
  for_each(flat_range(partition),
           [&partition, &part](int_t i) { partition[i] = int_t(part[i]); });

  return partition;
}

array<metis_idx_t, 1> compute_partitions_mesh(const Grid &grid, int n_parts) {
  auto n_cells = grid.n_cells;
  auto n_vertices = grid.n_vertices;
  auto max_neighbours = grid.max_neighbours;
  const auto &vertex_indices = grid.vertex_indices;

  auto ne = metis_idx_t(n_cells);
  auto nn = metis_idx_t(n_vertices);

  array<metis_idx_t, 1> eptr(n_cells + 1);
  array<metis_idx_t, 1> eind(shape_t<1>{n_cells * max_neighbours});

  for (int_t i = 0; i < n_cells; ++i) {
    eptr[i] = 3 * metis_idx_t(i);

    for (int_t k = 0; k < max_neighbours; ++k) {
      eind[integer_cast<int_t>(eptr[i]) + k]
          = integer_cast<metis_idx_t>(vertex_indices(i, k));
    }
  }
  eptr[n_cells] = 3 * metis_idx_t(n_cells);

  auto ncommon = metis_idx_t(max_neighbours - 1);
  auto nparts = metis_idx_t(n_parts);

  metis_idx_t objval;
  array<metis_idx_t, 1> epart(shape_t<1>{n_cells});
  array<metis_idx_t, 1> npart(shape_t<1>{n_vertices});

  // clang-format off
  METIS_PartMeshDual(&ne, &nn, eptr.raw(), eind.raw(), nullptr, nullptr,
                     &ncommon, &nparts, nullptr, nullptr, &objval, epart.raw(),
                     npart.raw());
  // clang-format on

  return epart;
}

array<int_t, 1>
compute_cell_permutation(const Grid &grid,
                         const array_const_view<int_t, 1> &cell_partition) {

  assert(grid.n_cells < int_t(std::numeric_limits<int>::max));

  auto n_cells = grid.n_cells;

  int n_mini_patches = zisa::max(2, int(n_cells) / 512);
  auto cell_mini_partitions = compute_partitions_mesh(grid, n_mini_patches);

  array<int_t, 1> sigma(n_cells);
  for_each(flat_range(sigma), [&sigma](int_t i) { sigma[i] = i; });

  std::sort(sigma.begin(),
            sigma.end(),
            [&cell_partition, &cell_mini_partitions](int_t i, int_t j) {
              if (cell_partition[i] == cell_partition[j]) {
                return cell_mini_partitions[i] < cell_mini_partitions[j];
              } else {
                return cell_partition[i] < cell_partition[j];
              }
            });

  return sigma;
}

array<int_t, 1>
compute_partition_boundaries(const array<int_t, 1> &cell_partition,
                             int n_parts) {

  array<int_t, 1> count(n_parts);
  fill(count, int_t(0));
  for_each(
      serial_policy{},
      PlainIndexRange(cell_partition.shape(0)),
      [&count, &cell_partition](int_t i) { count[cell_partition[i]] += 1; });

  array<int_t, 1> boundaries(n_parts + 1);
  boundaries(0) = 0;

  for (int p = 0; p < n_parts; ++p) {
    boundaries(p + 1) = boundaries(p) + count(p);
  }

  return boundaries;
}

array<int_t, 1> compute_inverse_permutation(const array<int_t, 1> &sigma) {
  auto tau = array<int_t, 1>(sigma.shape());
  for (int_t i = 0; i < sigma.shape(0); ++i) {
    tau(sigma(i)) = i;
  }

  return tau;
}

array<int_t, 2> renumbered_vertex_indices(const array<int_t, 2> &vertex_indices,
                                          const array<int_t, 1> &permutation) {
  auto vi = array<int_t, 2>(vertex_indices.shape());

  for (int_t i = 0; i < vertex_indices.shape(0); ++i) {
    for (int_t k = 0; k < vertex_indices.shape(1); ++k) {
      vi(i, k) = vertex_indices(permutation(i), k);
    }
  }

  return vi;
}

std::vector<std::vector<int_t>>
compute_effective_stencils(const array<StencilFamily, 1> &stencils) {
  std::vector<std::vector<int_t>> effective_stencils(stencils.size());
  for (int_t i = 0; i < stencils.size(); ++i) {
    effective_stencils[i] = stencils(i).local2global();
  }

  return effective_stencils;
}

PartitionedGrid compute_partitioned_grid(
    const Grid &grid, const array<StencilFamily, 1> &stencils, int_t n_parts) {

  std::vector<std::vector<int_t>> effective_stencils
      = compute_effective_stencils(stencils);

  auto cell_partition = compute_partition_full_stencil(
      grid, effective_stencils, integer_cast<int>(n_parts));

  auto sigma = compute_cell_permutation(grid, cell_partition);

  auto boundaries = compute_partition_boundaries(cell_partition,
                                                 integer_cast<int>(n_parts));

  return PartitionedGrid{
      std::move(cell_partition), std::move(boundaries), std::move(sigma)};
}

std::tuple<array<int_t, 2>, array<int_t, 1>>
extract_vertex_indices(const array<int_t, 2> &vertex_indices,
                       const std::function<bool(int_t)> &is_inside) {
  auto n_cells = vertex_indices.shape(0);
  auto max_vertices = vertex_indices.shape(1);

  int_t n_local_cells = 0;
  for (int_t i : index_range(n_cells)) {
    n_local_cells += (is_inside(i) ? 1 : 0);
  }

  auto local_vertex_indices
      = array<int_t, 2>(shape_t<2>{n_local_cells, vertex_indices.shape(1)});

  auto global_cell_indices = array<int_t, 1>(n_local_cells);

  int_t i_local = 0;
  for (int_t i : index_range(n_cells)) {
    if (is_inside(i)) {
      for (int_t k : index_range(max_vertices)) {
        local_vertex_indices(i_local, k) = vertex_indices(i, k);
      }
      global_cell_indices[i_local] = i;

      ++i_local;
    }
  }

  return {std::move(local_vertex_indices), std::move(global_cell_indices)};
}

StencilBasedIndicator::StencilBasedIndicator(
    const Grid &grid,
    const array<StencilFamily, 1> &stencils,
    const std::function<bool(int_t)> &is_interior) {
  LOG_ERR_IF(grid.n_cells != stencils.size(), "size mismatch.");

  for (int_t i : cell_indices(grid)) {
    if (is_interior(i)) {
      append(stencils[i].local2global());

      for (int_t k : neighbour_index_range(grid)) {
        if (grid.is_valid(i, k)) {
          int_t j = grid.neighbours(i, k);
          append(stencils[j].local2global());
        }
      }
    }
  }

  std::sort(good_indices.begin(), good_indices.end());
  auto end = std::unique(good_indices.begin(), good_indices.end());
  good_indices.resize(integer_cast<int_t>(end - good_indices.begin()));
}

bool StencilBasedIndicator::operator()(int_t i) const {
  return std::binary_search(good_indices.begin(), good_indices.end(), i);
}

void StencilBasedIndicator::append(const std::vector<int_t> &l2g) {
  good_indices.insert(good_indices.end(), l2g.begin(), l2g.end());
}

std::tuple<array<int_t, 2>, array<XYZ, 1>>
extract_vertices(const array<XYZ, 1> &vertices,
                 array<int_t, 2> local_vertex_indices) {

  std::vector<int_t> vi_local2old(local_vertex_indices.size());
  std::copy(local_vertex_indices.begin(),
            local_vertex_indices.end(),
            vi_local2old.begin());

  std::sort(vi_local2old.begin(), vi_local2old.end());

  auto n_unique = std::unique(vi_local2old.begin(), vi_local2old.end())
                  - vi_local2old.begin();
  vi_local2old.resize(integer_cast<size_t>(n_unique));

  std::map<int_t, int_t> vi_old2local;
  zisa::for_each(index_range(vi_local2old.size()),
                 [&](int_t i) { vi_old2local[vi_local2old[i]] = i; });

  zisa::for_each(index_range(local_vertex_indices.shape(0)), [&](int_t i) {
    for (int_t k : index_range(local_vertex_indices.shape(1))) {
      local_vertex_indices(i, k) = vi_old2local[local_vertex_indices(i, k)];
    }
  });

  int_t n_vertices_local = vi_old2local.size();
  array<XYZ, 1> local_vertices(n_vertices_local);

  zisa::for_each(index_range(local_vertices.shape(0)), [&](int_t i) {
    local_vertices(i) = vertices(vi_local2old[i]);
  });

  return {std::move(local_vertex_indices), std::move(local_vertices)};
}

std::tuple<array<int_t, 2>, array<XYZ, 1>, array<int_t, 1>>
extract_subgrid_v2(const Grid &grid,
                   const std::function<bool(int_t)> &is_inside) {
  return extract_subgrid_v2(grid.vertex_indices, grid.vertices, is_inside);
}

std::tuple<array<int_t, 2>, array<XYZ, 1>, array<int_t, 1>>
extract_subgrid_v2(const array<int_t, 2> &global_vertex_indices,
                   const array<XYZ, 1> &global_vertices,
                   const std::function<bool(int_t)> &is_inside) {

  auto [lvi, global_cell_indices]
      = extract_vertex_indices(global_vertex_indices, is_inside);

  auto [local_vertex_indices, local_vertices]
      = extract_vertices(global_vertices, std::move(lvi));

  return std::tuple{std::move(local_vertex_indices),
                    std::move(local_vertices),
                    std::move(global_cell_indices)};
}
DistributedGrid
extract_distributed_subgrid(const DistributedGrid &dgrid,
                            const array<int_t, 1> &global_indices) {

  return DistributedGrid{
      extract_subarray(dgrid.global_cell_indices, global_indices),
      extract_subarray(dgrid.partition, global_indices)};
}

std::map<int_t, int_t>
sparse_inverse_permutation(const array_const_view<int_t, 1> &sigma) {
  std::map<int_t, int_t> tau;
  for (int_t i = 0; i < sigma.size(); ++i) {
    tau[sigma[i]] = i;
  }
  tau[int_t(-1)] = int_t(-1);

  return tau;
}

array<StencilFamily, 1>
extract_stencils(const array<StencilFamily, 1> &global_stencils,
                 const array<int_t, 2> &global_neighbours,
                 const std::function<bool(int_t)> &global_is_inside,
                 const array<int_t, 1> &global_indices) {

  auto local_stencils = extract_subarray(global_stencils, global_indices);
  for (int_t i_local = 0; i_local < local_stencils.size(); ++i_local) {
    int_t i_global = global_indices[i_local];

    bool is_good = false;
    if (global_is_inside(i_global)) {
      is_good = true;
    } else {
      for (int_t k = 0; k < global_neighbours.shape(1); ++k) {
        int_t j_global = global_neighbours(i_global, k);

        if (j_global != int_t(-1) && global_is_inside(j_global)) {
          is_good = true;
          break;
        }
      }
    }

    if (!is_good) {
      local_stencils[i_local].truncate_to_first_order(i_global);
    }
  }

  auto tau = sparse_inverse_permutation(global_indices);
  for (int_t i = 0; i < local_stencils.size(); ++i) {
    local_stencils[i].apply_permutation(
        [&tau](int_t i_inner) { return tau.at(i_inner); });
  }

  return local_stencils;
}

// --- old extract subgrid -----------------------------------------------------
std::tuple<array<int_t, 2>, array<XYZ, 1>, array<int_t, 1>>
extract_subgrid(const Grid &grid,
                const PartitionedGrid &partitioned_grid,
                const array<StencilFamily, 1> &stencils,
                int_t k_part) {

  const auto &neighbours = grid.neighbours;
  const auto &sigma = partitioned_grid.permutation;
  auto tau = compute_inverse_permutation(sigma);
  const auto &boundaries = partitioned_grid.boundaries;
  const auto &partition = partitioned_grid.partition;

  auto n_cells = grid.n_cells;
  auto max_neighbours = grid.max_neighbours;
  auto n_cells_part = boundaries[k_part + 1] - boundaries[k_part];

  std::map<int_t, std::vector<int_t>> m;

  // Add the stencils of the interior.
  for (int_t i = 0; i < n_cells_part; ++i) {
    const auto &l2g = stencils[sigma(boundaries[k_part] + i)].local2global();

    for (auto j : l2g) {
      auto p = partition(j);
      if (p != k_part) {
        if (std::find(m[p].cbegin(), m[p].cend(), j) == m[p].cend()) {
          m[p].push_back(j);
        }
      }
    }
  }

  // Add stencils of the just-outside cells.
  for (int_t i = 0; i < n_cells_part; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {
      int_t j = neighbours(sigma(boundaries[k_part] + i), k);
      if (j < n_cells) {
        const auto &l2g = stencils[j].local2global();
        for (auto jj : l2g) {
          int_t p = partition(jj);
          if (p != k_part) {
            if (std::find(m[p].cbegin(), m[p].cend(), jj) == m[p].cend()) {
              m[p].push_back(jj);
            }
          }
        }
      }
    }
  }

  // count number of halo cells
  int_t n_halo = 0;
  for (const auto &[_, nbs] : m) {
    n_halo += nbs.size();
  }

  std::map<int_t, int_t> old2local;
  for (int_t i = 0; i < n_cells_part; ++i) {
    old2local[sigma(i + boundaries(k_part))] = i;
  }

  int_t i_local = n_cells_part;
  for (auto [p, halo] : m) {
    for (int_t i : halo) {
      old2local[i] = i_local;
      ++i_local;
    }
  }
  old2local[int_t(-1)] = int_t(-1);

  std::map<int_t, int_t> local2old;
  for (auto [old, local] : old2local) {
    local2old[local] = old;
  }

  int_t n_cells_halo = old2local.size() - 1 - n_cells_part;
  int_t n_cells_local = n_cells_part + n_cells_halo;

  array<int_t, 2> local_vertex_indices(
      {n_cells_part + n_cells_halo, max_neighbours});
  for (int_t i = 0; i < n_cells_local; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {
      local_vertex_indices(i, k) = grid.vertex_indices(local2old[i], k);
    }
  }

  std::vector<int_t> vi_local2old(local_vertex_indices.size());
  std::copy(local_vertex_indices.begin(),
            local_vertex_indices.end(),
            vi_local2old.begin());

  std::sort(vi_local2old.begin(), vi_local2old.end());

  auto n_unique = std::unique(vi_local2old.begin(), vi_local2old.end())
                  - vi_local2old.begin();
  vi_local2old.resize(integer_cast<size_t>(n_unique));

  std::map<int_t, int_t> vi_old2local
      = sparse_inverse_permutation(array_const_view(vi_local2old));

  for (auto &i : local_vertex_indices) {
    i = vi_old2local[i];
  }

  int_t n_vertices_local = vi_old2local.size();
  array<XYZ, 1> local_vertices(n_vertices_local);
  for (int_t i = 0; i < n_vertices_local; ++i) {
    local_vertices(i) = grid.vertices(vi_local2old[i]);
  }

  auto global_cell_indices = array<int_t, 1>(n_cells_local);
  for (int_t i = 0; i < n_cells_local; ++i) {
    global_cell_indices[i] = local2old[i];
  }

  return std::tuple{std::move(local_vertex_indices),
                    std::move(local_vertices),
                    std::move(global_cell_indices)};
}

}
