#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <tuple>

#if ZISA_HAS_METIS == 1
#include <metis.h>
#endif

#if ZISA_HAS_MPI == 1
#include <mpi.h>
#endif

#include <zisa/grid/grid.hpp>
#include <zisa/io/format_as_list.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/parallelization/domain_decomposition.hpp>
#include <zisa/reconstruction/stencil_family.hpp>

#include <thread>
#include <zisa/utils/timer.hpp>

namespace po = boost::program_options;

namespace zisa {
#if ZISA_HAS_METIS == 1
using metis_idx_t = ::idx_t;
#endif

std::mutex work_queue_mutex;

void save_partitioned_grid(const std::string &dirname,
                           const Grid &grid,
                           const StencilParams &stencil_params,
                           int n_parts,
                           int n_workers) {

  LOG_ERR_IF(n_parts <= 1, "You need 2 or more parts.");

  //  auto sharp_stencil_params = StencilFamilyParams{
  //      {3, 2, 2, 2, 2}, {"c", "b", "b", "b", "b"},
  //      {2.0, 1.5, 1.5, 1.5, 1.5}};
  //
  //  auto stencil_timer = Timer();
  //  auto effective_stencils
  //      = compute_effective_stencils(grid, sharp_stencil_params);
  //  PRINT(stencil_timer.elapsed_seconds());
  //
  //  auto metis_timer = Timer();
  //  //  auto partitioned_grid = compute_partitioned_grid(grid, n_parts);
  //  auto partitioned_grid
  //      = compute_partitioned_grid(grid, effective_stencils, n_parts);
  //  PRINT(metis_timer.elapsed_seconds());

  auto partitioned_grid = compute_partitioned_grid_by_sfc(grid, n_parts);

  const auto &permutation = partitioned_grid.permutation;
  const auto &partition = partitioned_grid.partition;
  const auto &boundaries = partitioned_grid.boundaries;

  std::vector<std::thread> workers;
  workers.reserve(n_workers);

  std::vector<int> work_queue(n_parts, 0);

  auto job = [&](int p) {
    auto [local_vertex_indices, local_vertices, global_cell_indices]
        = extract_subgrid(grid, partitioned_grid, stencil_params, p);

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
  };

  auto scheduling_loop = [&work_queue, job]() {
    while (true) {
      int p = -1;
      {
        auto lock = std::unique_lock(work_queue_mutex);

        int n_parts = work_queue.size();
        for (int pp = 0; pp < n_parts; ++pp) {
          if (work_queue[pp] == 0) {
            work_queue[pp] = 1;
            p = pp;
            break;
          }
        }

        if (p == -1) {
          return; // no more work to be had.
        }
      }

      job(p);
    }
  };

  for (int w = 0; w < n_workers; ++w) {
    workers.push_back(std::thread(scheduling_loop));
  }

  for (int w = 0; w < n_workers; ++w) {
    workers[w].join();
  }

  // auto filename = dirname + "/grid.h5";
  // auto hdf5_writer = HDF5SerialWriter(filename);
  // save(hdf5_writer, grid);
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
      ("workers,k", po::value<int>(), "Number of workers.")
      ;
  // clang-format on

  // first parse cmdline and check what config file to use
  po::store(po::parse_command_line(argc, argv, generic), options);

  if (options.count("help") != 0) {
    std::cout << generic << "\n";
    std::exit(EXIT_SUCCESS);
  }

  if (options.count("grid") == 0) {
    std::cout << "Missing argument `--gmsh GMSH`.\n";
    std::exit(EXIT_FAILURE);
  }

  if (options.count("partitions") == 0) {
    std::cout << "Missing argument `-n N`.\n";
    std::exit(EXIT_FAILURE);
  }

  if (options.count("partitions") == 0) {
    std::cout << "Missing argument `-k K`.\n";
    std::exit(EXIT_FAILURE);
  }

  if (options.count("output") == 0) {
    std::cout << "Missing argument `--output OUTPUT`.\n";
    std::exit(EXIT_FAILURE);
  }

  auto n_parts = options["partitions"].as<int>();
  auto n_workers = options["workers"].as<int>();
  auto gmsh_file = options["grid"].as<std::string>();
  auto part_file = options["output"].as<std::string>();
  auto grid = zisa::load_grid(gmsh_file);

  auto stencil_params = [&grid]() {
    if (grid->n_dims() == 2) {
      return zisa::StencilParams(5, "c", 8.0);
    } else if (grid->n_dims() == 3) {
      return zisa::StencilParams(4, "c", 8.0);
    }
    LOG_ERR("Broken logic.");
  }();

  zisa::save_partitioned_grid(
      part_file, *grid, stencil_params, n_parts, n_workers);

#if ZISA_HAS_MPI == 1
  MPI_Finalize();
#endif
}
