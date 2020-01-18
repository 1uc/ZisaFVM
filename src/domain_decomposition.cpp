#if ZISA_HAS_MPI == 1
#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <tuple>

#include <metis.h>
#include <mpi.h>

#include <zisa/grid/grid.hpp>
#include <zisa/io/format_as_list.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/loops/for_each.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/parallelization/halo_exchange.hpp>
#include <zisa/reconstruction/stencil_family.hpp>

#include <zisa/boundary/halo_exchange_bc.hpp>
#include <zisa/boundary/no_boundary_condition.hpp>

// I/O related includes
#include <zisa/io/dump_snapshot.hpp>
#include <zisa/io/gathered_visualization.hpp>
#include <zisa/parallelization/all_variables_gatherer.hpp>
#include <zisa/parallelization/mpi_single_node_array_gatherer.hpp>

#include <zisa/parallelization/domain_decomposition.hpp>

namespace po = boost::program_options;

namespace zisa {
using metis_idx_t = ::idx_t;

template <class Int = int_t>
class TopologicalGridView {
public:
  array_const_view<Int, 2> vertex_indices;
  array_const_view<Int, 2> neighbours;
  Int n_cells;
  Int n_vertices;
  Int max_neighbours;

public:
  TopologicalGridView(const array_const_view<Int, 2> &vertex_indices,
                      const array_const_view<Int, 2> &neighbours,
                      Int n_vertices,
                      Int max_neighbours)
      : vertex_indices(vertex_indices),
        neighbours(neighbours),
        n_cells(vertex_indices.shape(0)),
        n_vertices(n_vertices),
        max_neighbours(max_neighbours) {}
};

TopologicalGridView<int_t> make_topological_grid_view(const Grid &grid) {
  return {grid.vertex_indices,
          grid.neighbours,
          grid.n_vertices,
          grid.max_neighbours};
}

void save_partitioned_grid(const std::string &filename,
                           const Grid &grid,
                           const array<StencilFamily, 1> &stencils,
                           int n_parts) {

  auto partitioned_grid = compute_partitioned_grid(grid, stencils, n_parts);
  const auto &permutation = partitioned_grid.permutation;

  auto vertex_indices
      = renumbered_vertex_indices(grid.vertex_indices, permutation);

  auto renumbered_grid
      = Grid(GMSHElementType::triangle, grid.vertices, vertex_indices, 1);

  const auto &partition = partitioned_grid.partition;
  const auto &boundaries = partitioned_grid.boundaries;

  int_t mpi_rank = zisa::mpi::rank(MPI_COMM_WORLD);
  {
    auto [local_vertex_indices, local_vertices, local_stencils, halo]
        = extract_subgrid(grid, partitioned_grid, stencils, mpi_rank);

    auto element_type = grid.max_neighbours == 3
                            ? zisa::GMSHElementType::triangle
                            : zisa::GMSHElementType::tetrahedron;
    int_t quad_deg = 1;

    auto local_grid = zisa::Grid(
        element_type, local_vertices, local_vertex_indices, quad_deg);

    auto n_cells_local = local_grid.n_cells;
    int_t n_owned_cells = partitioned_grid.boundaries[mpi_rank + 1]
                          - partitioned_grid.boundaries[mpi_rank];

    auto data_local = GridVariables({n_cells_local, 5});
    zisa::fill(data_local, -123.0);
    for(int_t i = 0; i < n_owned_cells; ++i) {
      for(int_t k = 0; k < data_local.shape(1); ++k) {
        data_local(i, k) = double(mpi_rank);
      }
    }

    auto array_info
        = std::make_shared<DistributedArrayInfo>(partitioned_grid.boundaries);
    auto cvars_gatherer
        = std::make_unique<MPISingleNodeArrayGatherer<double, 2>>(
            array_info, MPI_COMM_WORLD, 388932);
    auto avars_gatherer
        = std::make_unique<MPISingleNodeArrayGatherer<double, 2>>(
            array_info, MPI_COMM_WORLD, 333292);
    auto all_vars_gatherer = std::make_shared<AllVariablesGatherer>(
        std::move(cvars_gatherer), std::move(avars_gatherer), n_owned_cells);

    auto all_vars_local = AllVariables({n_cells_local, 5, 0});
    all_vars_local.cvars = data_local;

    auto halo_exchange = std::make_shared<MPIHaloExchange>(
        make_mpi_halo_exchange(halo, MPI_COMM_WORLD));

    auto nobc = std::make_shared<NoBoundaryCondition>();
    auto bc = std::make_shared<HaloExchangeBC>(nobc, halo_exchange);
    bc->apply(all_vars_local, 0.0);

    auto euler = std::make_shared<Euler<IdealGasEOS, NoGravity>>();
    auto fng = std::make_shared<FileNameGenerator>("f", "%d", ".h5");

    std::shared_ptr<Visualization> euler_visualization = nullptr;
    if (mpi_rank == 0) {
      euler_visualization
          = std::make_shared<SerialDumpSnapshot<Euler<IdealGasEOS, NoGravity>>>(
              euler, fng);
    }

    auto gathered_visualization = std::make_shared<GatheredVisualization>(
        std::move(all_vars_gatherer),
        std::move(euler_visualization),
        AllVariablesDimensions{grid.n_cells, 5ul, 0ul});

    auto simulation_clock = SerialSimulationClock(nullptr, nullptr);

    (*gathered_visualization)(all_vars_local, simulation_clock);

    if (mpi_rank == 0) {
      auto writer = zisa::HDF5SerialWriter("grid_full.h5");
      save(writer, renumbered_grid);
    }

    {
      auto writer = zisa::HDF5SerialWriter(
          string_format("grid_part-%04d.h5", mpi_rank));
      save(writer, local_grid);
    }

    {
      auto writer = zisa::HDF5SerialWriter(
          string_format("data_part-%04d.h5", mpi_rank));
      save(writer, all_vars_local.cvars, "upsilon");
    }
  }
}
}

int main(int argc, char *argv[]) {
  int requested = MPI_THREAD_MULTIPLE;
  int provided = -1;

  MPI_Init_thread(&argc, &argv, requested, &provided);
  LOG_ERR_IF(requested > provided,
             "MPI does not support the requested level of multi-threading.");

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

  auto stencil_params = zisa::StencilFamilyParams(
      {3, 2, 2, 2}, {"c", "b", "b", "b"}, {5.0, 2.0, 2.0, 2.0});

  auto stencils = zisa::compute_stencil_families(*grid, stencil_params);

  zisa::save_partitioned_grid(part_file, *grid, stencils, n_parts);

  MPI_Finalize();
}

#else
#include <zisa/config.hpp>
int main() {
  LOG_ERR("Must be compiled with `-DZISA_HAS_MPI=1`.");
}
#endif
