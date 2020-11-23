#include <zisa/io/file_name_generator.hpp>
#include <zisa/io/load_snapshot.hpp>
#include <zisa/io/parallel_dump_snapshot.hpp>
#include <zisa/io/parallel_load_snapshot.hpp>
#include <zisa/model/ideal_gas_eos.hpp>
#include <zisa/model/local_eos_state.hpp>
#include <zisa/mpi/io/gathered_vis_info.hpp>
#include <zisa/mpi/io/gathered_visualization_factory.hpp>
#include <zisa/mpi/io/scattered_data_source_factory.hpp>
#include <zisa/mpi/math/distributed_reference_solution.hpp>
#include <zisa/mpi/parallelization/mpi_halo_exchange.hpp>
#include <zisa/mpi/parallelization/mpi_single_node_array_gatherer_decl.hpp>
#include <zisa/mpi/parallelization/mpi_single_node_array_scatterer_decl.hpp>
#include <zisa/ode/simulation_clock.hpp>
#include <zisa/parallelization/local_grid.hpp>

namespace zisa {

void test_distributed_reference() {
  const auto &world_comm = MPI_COMM_WORLD;
  auto world_rank = zisa::mpi::rank(world_comm);
  auto world_comm_size = zisa::mpi::size(world_comm);
  LOG_ERR_IF(world_comm_size != 4, "This test needs exactly 4 MPI tasks.");

  int_t n_dummy = 0;
  int_t n_cvars = 5;
  int_t n_avars = 0;

  auto simulation_clock = std::make_shared<SerialSimulationClock>(
      std::make_shared<DummyTimeKeeper>(),
      std::make_shared<DummyPlottingSteps>());

  auto all_var_dims = AllVariablesDimensions{n_dummy, n_cvars, n_avars};

  auto serialize = [all_var_dims,
                    simulation_clock](const std::string &filename,
                                      MPI_Comm small_comm,
                                      const DistributedGrid &small_dgrid,
                                      const AllVariables &all_vars) {
    int n_writers = 2;
    auto vis_info = make_gathered_vis_info(small_comm, small_dgrid, n_writers);
    auto gatherer_factory
        = make_mpi_single_node_array_gatherer_factory(*vis_info);

    auto fng = std::make_shared<SingleFileNameGenerator>(filename);
    auto file_dims = make_hdf5_unstructured_file_dimensions(*vis_info);

    auto local_eos
        = std::make_shared<LocalEOSState<IdealGasEOS>>(/* gamma = */ 1.2,
                                                       /* R = */ 1.0);

    auto dump_snapshot = std::make_shared<ParallelDumpSnapshot<IdealGasEOS>>(
        local_eos, fng, file_dims);

    auto vis = make_gathered_visualization(std::move(vis_info),
                                           std::move(gatherer_factory),
                                           std::move(dump_snapshot),
                                           all_var_dims);

    (*vis)(all_vars, *simulation_clock);
  };

  auto subgrid_name = string_format(
      "grids/unit_tests/reference_soln/partitioned/4/subgrid-%04d", world_rank);

  auto stencil_params = StencilFamilyParams({});
  auto qr_degrees = QRDegrees{2, 2, 2};

  auto large_grid = std::get<2>(
      load_local_grid(subgrid_name, stencil_params, qr_degrees, world_rank));

  auto interpolation = [](int_t, const XYZ &x, int_t) {
    return sin(2.0 * zisa::pi * zisa::norm(x));
  };

  auto n_vars = n_cvars + n_avars;
  auto dref = DistributedReferenceSolution(
      serialize, large_grid, interpolation, n_vars);
}

AllVariables
reference_all_vars(const Grid &grid, int_t n_cvars, int_t n_avars) {
  auto n_cells = grid.n_cells;

  auto all_vars = AllVariables({n_cells, n_cvars, n_avars});
  for (int_t i = 0; i < n_cells; ++i) {
    auto x = zisa::norm(grid.cell_centers[i]);

    for (int_t k = 0; k < n_cvars; ++k) {
      all_vars.cvars(i, k) = 10.0 * (k + 1) + zisa::sin2pi(x);
    }

    for (int_t k = 0; k < n_avars; ++k) {
      all_vars.avars(i, k) = 10.0 * (n_cvars + k + 1) + zisa::sin2pi(x);
    }
  }

  return all_vars;
}

void test_distributed_write() {
  // Stage 0: setup stuff.
  const auto &world_comm = MPI_COMM_WORLD;
  auto world_rank = zisa::mpi::rank(world_comm);
  auto world_size = zisa::mpi::size(world_comm);

  int n_writers = 2;
  int_t n_dummy = 0;
  int_t n_cvars = 5;
  int_t n_avars = 2;
  double t_final = 12.0;
  int_t n_steps = 14;
  double atol = 1e-10;

  double gamma = 1.2;
  double rgas = 1.0;

  LOG_ERR_IF(world_size != 4, "This must be run with 4 MPI tasks.");

  auto filename = std::string("distributed_data.h5");
  auto grid_directory = std::string("grids/unit_tests/distributed_io");
  auto subgrid_name = string_format("%s/partitioned/4/subgrid-%04d.msh.h5",
                                    grid_directory.c_str(),
                                    world_rank);

  auto stencil_params = StencilFamilyParams{
      {3, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}};

  auto qr_degrees = QRDegrees{2, 2, 4};

  auto [stencils, dgrid, grid]
      = load_local_grid(subgrid_name, stencil_params, qr_degrees, world_rank);

  auto vis_info = make_gathered_vis_info(world_comm, *dgrid, n_writers);
  auto gatherer_factory
      = make_mpi_single_node_array_gatherer_factory(*vis_info);

  auto all_var_dims = AllVariablesDimensions{n_dummy, n_cvars, n_avars};
  auto fng = std::make_shared<SingleFileNameGenerator>(filename);
  auto file_dims = make_hdf5_unstructured_file_dimensions(*vis_info);

  // Stage 1: write data to disk.
  auto local_eos = std::make_shared<LocalEOSState<IdealGasEOS>>(gamma, rgas);
  auto dump_snapshot = std::make_shared<ParallelDumpSnapshot<IdealGasEOS>>(
      local_eos, fng, file_dims);

  auto out_clock = std::make_shared<SerialSimulationClock>(
      std::make_shared<DummyTimeKeeper>(),
      std::make_shared<DummyPlottingSteps>());
  out_clock->advance_to(t_final, n_steps);

  auto visualization = make_gathered_visualization(
      vis_info, gatherer_factory, dump_snapshot, all_var_dims);

  auto out_vars = reference_all_vars(*grid, n_cvars, n_avars);
  (*visualization)(out_vars, *out_clock);

  visualization->wait();
  MPI_Barrier(world_comm);

  // Stage 2: load the array back into memory.
  auto scatterer_factory
      = make_mpi_single_node_array_scatterer_factory(*vis_info);

  {
    int_t n_cells = grid->n_cells;
    auto in_dims = AllVariablesDimensions{n_cells, n_cvars, n_avars};

    auto in_vars = AllVariables(in_dims);
    auto in_clock = std::make_shared<SerialSimulationClock>(
        std::make_shared<DummyTimeKeeper>(),
        std::make_shared<DummyPlottingSteps>());

    auto load_snapshot = std::make_shared<ParallelLoadSnapshot>(fng, file_dims);
    auto halo_exchange = std::make_shared<MPIHaloExchange>(
        make_mpi_halo_exchange(*dgrid, world_comm));

    auto data_source = make_scattered_data_source(vis_info,
                                                  scatterer_factory,
                                                  load_snapshot,
                                                  halo_exchange,
                                                  all_var_dims);

    (*data_source)(in_vars, *in_clock);

    LOG_ERR_IF(in_vars.dims() != in_dims, "Wrong dimensions.");

    for (int_t i = 0; i < n_cells; ++i) {
      for (int_t k = 0; k < n_cvars; ++k) {
        auto is_good
            = zisa::abs(in_vars.cvars(i, k) - out_vars.cvars(i, k)) < atol;

        PRINT_IF(!is_good, in_vars.cvars(i, k) - out_vars.cvars(i, k));
        PRINT_IF(!is_good, world_rank);
        LOG_ERR_IF(!is_good, "Failed.");
      }

      for (int_t k = 0; k < n_avars; ++k) {
        auto is_good
            = zisa::abs(in_vars.avars(i, k) - out_vars.avars(i, k)) < atol;

        PRINT_IF(!is_good, in_vars.avars(i, k) - out_vars.avars(i, k));
        PRINT_IF(!is_good, world_rank);
        LOG_ERR_IF(!is_good, "Failed.");
      }
    }
  }

  PRINT("SUCCESS!");
}

}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  zisa::test_distributed_write();

  MPI_Finalize();
}
