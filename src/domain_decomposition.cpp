#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <tuple>

#include <metis.h>

#if ZISA_HAS_MPI == 1
#include <mpi.h>
#endif

#include <zisa/grid/grid.hpp>
#include <zisa/io/format_as_list.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/parallelization/domain_decomposition.hpp>
#include <zisa/reconstruction/stencil_family.hpp>

namespace po = boost::program_options;

namespace zisa {
using metis_idx_t = ::idx_t;

void save_partitioned_grid(const std::string &dirname,
                           const Grid &grid,
                           const array<StencilFamily, 1> &stencils,
                           int n_parts) {

  auto partitioned_grid = compute_partitioned_grid(grid, stencils, n_parts);
  const auto &permutation = partitioned_grid.permutation;

  const auto &partition = partitioned_grid.partition;
  const auto &boundaries = partitioned_grid.boundaries;

  for (int p = 0; p < n_parts; ++p) {
    auto [local_vertex_indices, local_vertices, global_cell_indices]
        = extract_subgrid(grid, partitioned_grid, stencils, p);

    int_t n_cells_local = local_vertex_indices.shape(0);

    auto local_partition = array<int_t, 1>(n_cells_local);
    for (int_t i = 0; i < n_cells_local; ++i) {
      local_partition[i] = partition[global_cell_indices[i]];
    }

    auto filename = dirname + string_format("/subgrid-%04d.msh.h5", p);
    auto hdf5_writer = HDF5SerialWriter(filename);

    hdf5_writer.write_scalar(grid.n_dims(), "n_dims");
    save(hdf5_writer, local_vertex_indices, "vertex_indices");
    save(hdf5_writer, local_vertices, "vertices");
    save(hdf5_writer, local_partition, "partition");
    save(hdf5_writer, global_cell_indices, "global_cell_indices");
  }

  auto filename = dirname + "/grid.h5";
  auto hdf5_writer = HDF5SerialWriter(filename);
  save(hdf5_writer, grid);
}
}

int main(int argc, char *argv[]) {
#if ZISA_HAS_MPI == 1
  MPI_Init(&argc, &argv);
#endif

  po::variables_map options;

  // generic options
  po::options_description generic("Generic options");

  // clang-format off
  generic.add_options()
      ("help,h", "produce this message")
      ("grid", po::value<std::string>(), "GMSH grid to partition.")
      ("output,o", po::value<std::string>(), "Name of the folder in which to place the partitioned grids.")
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
  auto grid = zisa::load_grid(gmsh_file);

  auto stencil_params = [&grid]() {
    if (grid->n_dims() == 2) {
      return zisa::StencilFamilyParams({5}, {"c"}, {8.0});
    } else if (grid->n_dims() == 3) {
      return zisa::StencilFamilyParams({4}, {"c"}, {8.0});
    }
    LOG_ERR("Broken logic.");
  }();

  auto stencils = zisa::compute_stencil_families(*grid, stencil_params);

  zisa::save_partitioned_grid(part_file, *grid, stencils, n_parts);

#if ZISA_HAS_MPI == 1
  MPI_Finalize();
#endif
}
