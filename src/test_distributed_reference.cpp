#include <zisa/io/file_name_generator.hpp>
#include <zisa/io/gathered_vis_info.hpp>
#include <zisa/io/parallel_dump_snapshot.hpp>
#include <zisa/model/ideal_gas_eos.hpp>
#include <zisa/model/local_eos_state.hpp>
#include <zisa/mpi/io/gathered_visualization_factory.hpp>
#include <zisa/mpi/math/distributed_reference_solution.hpp>
#include <zisa/mpi/parallelization/mpi_single_node_array_gatherer_decl.hpp>
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

}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  zisa::test_distributed_reference();

  MPI_Finalize();
}
