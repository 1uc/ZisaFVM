#include <boost/program_options.hpp>
#include <iostream>
#include <string>

#include <metis.h>

#include <zisa/grid/grid.hpp>
#include <zisa/io/format_as_list.hpp>
#include <zisa/loops/for_each.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>

namespace po = boost::program_options;

namespace zisa {
using metis_idx_t = ::idx_t;

template <class Int = int_t>
class TopologicalGridView {
public:
  array_view<Int, 2> vertex_indices;
  Int n_cells;
  Int n_vertices;
  Int max_neighbours;

public:
  TopologicalGridView(const array_view<Int, 2> &vertex_indices,
                      Int n_vertices,
                      Int max_neighbours)
      : vertex_indices(vertex_indices),
        n_cells(vertex_indices.shape(0)),
        n_vertices(n_vertices),
        max_neighbours(max_neighbours) {}
};

array<metis_idx_t, 1> compute_partitions(const TopologicalGridView<int_t> &grid,
                                         int n_parts) {
  auto n_cells = grid.n_cells;
  auto n_vertices = grid.n_vertices;
  auto max_neighbours = grid.max_neighbours;
  const auto &vertex_indices = grid.vertex_indices;

  auto ne = metis_idx_t(n_cells);
  auto nn = metis_idx_t(n_vertices);

  array<metis_idx_t, 1> eptr(n_cells + 1);
  array<metis_idx_t, 1> eind({n_cells * max_neighbours});

  for (int_t i = 0; i < n_cells; ++i) {
    eptr[i] = 3 * metis_idx_t(i);

    for (int_t k = 0; k < max_neighbours; ++k) {
      eind[eptr[i] + k] = vertex_indices(i, k);
    }
  }
  eptr[n_cells] = 3 * metis_idx_t(n_cells);

  auto ncommon = metis_idx_t(max_neighbours);
  auto nparts = metis_idx_t(n_parts);

  metis_idx_t objval;
  array<metis_idx_t, 1> epart({n_cells});
  array<metis_idx_t, 1> npart({n_vertices});

  // clang-format off
  METIS_PartMeshDual(&ne, &nn, eptr.raw(), eind.raw(), nullptr, nullptr,
                     &ncommon, &nparts, nullptr, nullptr, &objval, epart.raw(),
                     npart.raw());
  // clang-format on

  return epart;
}

array<int_t, 1>
compute_cell_permutation(const TopologicalGridView<int_t> &grid,
                         const array<metis_idx_t, 1> &cell_partition) {

  assert(grid.n_cells < int_t(std::numeric_limits<int>::max));

  auto n_cells = grid.n_cells;
  auto cell_mini_partitions = compute_partitions(grid, int(n_cells) / 16);

  array<int_t, 1> sigma(n_cells);

  for_each(PlainIndexRange(0, grid.n_cells),
           [&sigma](int_t i) { sigma[i] = i; });

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
compute_partition_boundaries(const array<metis_idx_t, 1> &cell_partition,
                             int n_parts) {

  array<int_t, 1> count(n_parts);
  fill(count, int_t(0));
  for_each(
      serial_policy{},
      PlainIndexRange(cell_partition.shape(0)),
      [&count, &cell_partition](int_t i) { count[cell_partition[i]] += 1; });

  array<int_t, 1> boundaries(n_parts + 1);
  boundaries(0) = 0;

  PRINT(format_as_list(count));

  for (int p = 0; p < n_parts; ++p) {
    boundaries(p+1) = boundaries(p) + count(p);
  }

  return boundaries;
}

std::pair<array<int_t, 2>, array<int_t, 1>>
compute_partitioned_vertex_indices(const TopologicalGridView<int_t> &grid,
                                   int n_parts) {

  auto cell_partition = compute_partitions(grid, n_parts);
  auto sigma = compute_cell_permutation(grid, cell_partition);

  auto n_cells = grid.n_cells;
  auto max_vertices = grid.max_neighbours;

  PRINT(n_cells);

  array<int_t, 2> pvi({n_cells, max_vertices});

  for_each(PlainIndexRange(n_cells), [&pvi, &sigma, &grid](int_t i) {
    const auto &vi = grid.vertex_indices;
    auto max_vertices = grid.max_neighbours;
    for (int_t k = 0; k < max_vertices; ++k) {
      pvi(i, k) = vi(sigma(i), k);
    }
  });

  auto boundaries = compute_partition_boundaries(cell_partition, n_parts);
  PRINT(format_as_list(boundaries));

  return {std::move(pvi), std::move(boundaries)};
}

void save_partitioned_grid(const std::string &filename,
                           Grid &unpartitioned_grid,
                           int n_parts) {

  auto grid_view = TopologicalGridView<int_t>(
      unpartitioned_grid.vertex_indices,
      unpartitioned_grid.n_vertices,
      unpartitioned_grid.max_neighbours);

  auto [vertex_indices, boundaries]
      = compute_partitioned_vertex_indices(grid_view, n_parts);

  auto writer = HDF5SerialWriter(filename);
  save(writer, vertex_indices, "vertex_indices");
  save(writer, unpartitioned_grid.vertices, "vertices");
  writer.open_group("domain_decomposition");
  save(writer, boundaries, "boundaries");
}

}

int main(int argc, char *argv[]) {
  po::variables_map options;

  // generic options
  po::options_description generic("Generic options");

  // clang-format off
  generic.add_options()
      ("help,h", "produce this message")
      ("grid", po::value<std::string>(), "GMSH grid to partition.")
      ("output,o", po::value<std::string>(), "Name of the partitioned grid.")
      ("partitions,n", po::value<int>(), "Number of partitions.")
      ;
  // clang-format on

  // first parse cmdline and check what config file to use
  po::store(po::parse_command_line(argc, argv, generic), options);

  if (options.count("help") != 0) {
    std::cout << generic << "\n";
    std::exit(EXIT_SUCCESS);
  }

  if (options.count("grid") == 0) {
    std::cout << "Missing argument `--grid GRID`.\n";
    std::exit(EXIT_FAILURE);
  }

  if (options.count("partitions") == 0) {
    std::cout << "Missing argument `-n N`.\n";
    std::exit(EXIT_FAILURE);
  }

  if (options.count("output") == 0) {
    std::cout << "Missing argument `--output OUTPUT`.\n";
    std::exit(EXIT_FAILURE);
  }

  auto n_parts = options["partitions"].as<int>();
  auto gmsh_file = options["grid"].as<std::string>();
  auto part_file = options["output"].as<std::string>();
  auto grid = zisa::load_gmsh(gmsh_file);

  auto grid_writer = zisa::HDF5SerialWriter("unpartitioned_grid.h5");
  zisa::save(grid_writer, *grid);

  zisa::save_partitioned_grid(part_file, *grid, n_parts);
}
