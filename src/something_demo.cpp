#if ZISA_HAS_MPI == 1
#include <iostream>
#include <random>
#include <string>
#include <tuple>

#include <zisa/datetime.hpp>
#include <zisa/mpi/io/hdf5_unstructured_writer.hpp>
#include <zisa/mpi/mpi.hpp>

#include <boost/program_options.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/io/gathered_visualization.hpp>
#include <zisa/io/parallel_dump_snapshot.hpp>
#include <zisa/model/euler_factory.hpp>
#include <zisa/mpi/parallelization/mpi_single_node_array_gatherer.hpp>
#include <zisa/parallelization/all_variables_gatherer.hpp>
#include <zisa/parallelization/distributed_grid.hpp>

namespace po = boost::program_options;

namespace zisa {
void test_parallel_hdf5(const std::string &output_file) {
  auto mpi_rank = zisa::mpi::rank();
  auto comm_size = zisa::mpi::size();
  size_t n_total_cells = (239000 / comm_size) * comm_size;

  LOG_ERR_IF(n_total_cells % comm_size != 0, "incorrect number of processes.");

  size_t n_local_cells = n_total_cells / comm_size;
  auto gids = std::vector<hsize_t>(n_total_cells);
  for (size_t i = 0; i < n_total_cells; ++i) {
    gids[i] = i;
  }

  std::swap(gids[8], gids[n_total_cells - 8]);
  std::swap(gids[n_local_cells + 7], gids[n_total_cells - n_local_cells - 9]);

  //  std::mt19937 md((unsigned int)(1298375861));
  //  std::shuffle(gids.begin(), gids.end(), md);

  auto ids = std::vector<hsize_t>(n_local_cells);
  std::copy(gids.begin() + mpi_rank * n_local_cells,
            gids.begin() + (mpi_rank + 1) * n_local_cells,
            ids.begin());

  auto upsilon = array<double, 1>(ids.size());
  for (auto i : index_range(ids.size())) {
    upsilon(i) = ids[i];
  }

  auto tau = array<double, 2>(shape_t<2>{ids.size() + 42, 3});
  for (auto i : index_range(ids.size())) {
    tau(i, 0) = 100 * ids[i];
    tau(i, 1) = 100 * ids[i] + 1;
    tau(i, 2) = 100 * ids[i] + 2;
  }

  auto file_dims
      = make_hdf5_unstructured_file_dimensions(0, ids, MPI_COMM_WORLD);

  zisa::mpi::barrier(MPI_COMM_WORLD);
  auto t_start = zisa::current_time_stamp();
  {
    auto writer = HDF5UnstructuredWriter(output_file, file_dims);
    save(writer, upsilon, "upsilon");
    //  save(writer, tau, "tau");
  }
  auto elapsed = zisa::elapsed_seconds_since(t_start);
  PRINT(string_format("%d t_elapsed = %f s", mpi_rank, elapsed));
}

void test_gathered_hdf5() {
  auto world_comm = MPI_COMM_WORLD;

  auto world_size = zisa::mpi::size();
  auto world_rank = zisa::mpi::rank();

  auto grid_file = string_format(
      "grids/gaussian_bump-1/%d/subgrid-%04d.msh.h5", world_size, world_rank);

  auto grid = zisa::load_grid(grid_file);
  auto dgrid = zisa::load_distributed_grid(grid_file);

  auto n_local_cells = int(std::count(
      dgrid.partition.begin(), dgrid.partition.end(), int_t(world_rank)));

  int n_total_cells = -1;

  MPI_Allreduce(
      &n_local_cells, &n_total_cells, 1, MPI_INT, MPI_SUM, world_comm);

  PRINT(n_total_cells);

  auto ideal_vis_comm_size = 4;
  auto pes_per_block
      = (world_size + ideal_vis_comm_size - 1) / ideal_vis_comm_size;

  auto part_id = world_rank / pes_per_block;
  auto vis_comm = zisa::mpi::comm_split(world_comm, part_id, world_rank);
  auto vis_size = zisa::mpi::size(vis_comm);
  auto vis_rank = zisa::mpi::rank(vis_comm);
  auto vis_tag = 1948293;

  auto h5_comm = zisa::mpi::comm_split(world_comm, vis_rank, world_rank);

  std::shared_ptr<GatheredVisualization> gathered_visualization;

  if (vis_rank == 0) {
    PRINT(string_format("%d %d %d\n", world_rank, part_id, vis_rank));

    auto vis_cells_per_pe = array<int_t, 1>(shape_t<1>{vis_size});
    vis_cells_per_pe[0] = n_local_cells;
    zisa::mpi::gather(array_view(vis_cells_per_pe), 0, vis_comm);
    auto n_vis_cells = std::reduce(vis_cells_per_pe.begin(),
                                   vis_cells_per_pe.end(),
                                   int_t(0),
                                   [](int_t i, int_t j) { return i + j; });

    auto gids = array<int_t, 1>(shape_t<1>{n_vis_cells});
    auto vis_boundaries = array<int_t, 1>(shape_t<1>{vis_size + 1});
    vis_boundaries[0] = 0;
    for (int_t i = 0; i < vis_size; ++i) {
      vis_boundaries[i + 1] = vis_boundaries[i] + vis_cells_per_pe[i];
    }

    for (int_t i = 0; i < n_local_cells; ++i) {
      gids[i] = dgrid.global_cell_indices[i];
    }

    for (int_t p = 1; p < vis_size; ++p) {
      auto i0 = vis_boundaries[p];
      auto i1 = vis_boundaries[p + 1];
      zisa::mpi::recv(slice(array_view(gids), i0, i1), p, vis_tag, vis_comm);
    }

    //    auto sigma = std::make_shared<Permutation>(
    //        Permutation{array<Cycle, 1>(shape_t<1>{0})});
    auto sigma_array = array<int_t, 1>(gids.shape());
    for (int_t i = 0; i < gids.shape(0); ++i) {
      sigma_array[i] = i;
    }
    std::sort(sigma_array.begin(),
              sigma_array.end(),
              [&gids](size_t i, size_t j) { return gids[i] < gids[j]; });

    auto sigma = std::make_shared<Permutation>(
        factor_permutation(array_const_view(sigma_array)));
    apply_permutation(array_view(gids), *sigma);

    std::vector<hsize_t> hids(gids.size());
    for (int_t i = 0; i < hids.size(); ++i) {
      if (i > 0) {
        LOG_ERR_IF(gids[i - 1] >= gids[i], "Failed monotonicity.");
      }
      hids[i] = integer_cast<hsize_t>(gids[i]);
    }

    auto file_dims = make_hdf5_unstructured_file_dimensions(0, hids, h5_comm);

    auto darray_info = make_distributed_array_info(vis_boundaries);

    auto gatherer_factory
        = MPISingleNodeArrayGathererFactory(darray_info, vis_comm, vis_tag + 1);

    // Careful, order in which the gatherers are created matters.
    // Therefore, we can't create inline in the AllVariablesGatherer
    // constructor.
    auto cvars_gatherer = gatherer_factory.create_pointer<double, 2>();
    auto avars_gatherer = gatherer_factory.create_pointer<double, 2>();
    auto all_var_gatherer
        = std::make_unique<AllVariablesGatherer>(std::move(cvars_gatherer),
                                                 std::move(avars_gatherer),
                                                 n_local_cells,
                                                 n_vis_cells);

    using euler_t = Euler<IdealGasEOS, ConstantGravityRadial>;
    auto euler = std::make_shared<euler_t>(make_default_euler());

    auto fng = make_file_name_generator("", "dio-", "%04d", ".h5");

    auto dump_snapshot = std::make_shared<ParallelDumpSnapshot<euler_t>>(
        euler, fng, file_dims);

    auto all_var_dims = AllVariablesDimensions{n_vis_cells, 5, 0};

    gathered_visualization = std::make_shared<GatheredVisualization>(
        std::move(all_var_gatherer), sigma, dump_snapshot, all_var_dims);

  } else {
    int_t n_local_cells_ = n_local_cells;
    zisa::mpi::gather(array_view(shape_t<1>{1}, &n_local_cells_), 0, vis_comm);
    zisa::mpi::send(const_slice(array_const_view(dgrid.global_cell_indices),
                                0,
                                n_local_cells),
                    0,
                    vis_tag,
                    vis_comm);

    auto gatherer_factory
        = MPISingleNodeArrayGathererFactory(nullptr, vis_comm, vis_tag + 1);

    auto cvars_gatherer = gatherer_factory.create_pointer<double, 2>();
    auto avars_gatherer = gatherer_factory.create_pointer<double, 2>();
    auto all_var_gatherer
        = std::make_unique<AllVariablesGatherer>(std::move(cvars_gatherer),
                                                 std::move(avars_gatherer),
                                                 n_local_cells,
                                                 int_t(-1));

    gathered_visualization = std::make_shared<GatheredVisualization>(
        std::move(all_var_gatherer),
        nullptr,
        nullptr,
        AllVariablesDimensions{0, 0, 0});
  }

  AllVariablesDimensions dims{grid->n_cells, 5, 0};
  auto u0 = AllVariables(dims);

  for (int_t i = 0; i < grid->n_cells; ++i) {
    u0.cvars(i, 0) = world_rank;
    u0.cvars(i, 1) = zisa::sin2pi(zisa::norm(grid->cell_centers[i]));
    u0.cvars(i, 2) = 0.0;
    u0.cvars(i, 3) = 0.0;
    u0.cvars(i, 4) = dgrid.global_cell_indices[i];
  }

  auto sim_clock = SerialSimulationClock(nullptr, nullptr);
  zisa::mpi::barrier(world_comm);
  auto t_start = zisa::current_time_stamp();
  (*gathered_visualization)(u0, sim_clock);
  (*gathered_visualization).wait();
  auto elapsed = zisa::elapsed_seconds_since(t_start);
  PRINT(string_format("%d t_elapsed = %f s", world_rank, elapsed));
}

void test_open_close_open() {
  auto filename = std::string("foo.h5");
  auto h5_file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  H5Fclose(h5_file);

  h5_file
      = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  H5Fclose(h5_file);

  { auto reader = HDF5SerialReader(filename); }
  { auto writer = HDF5SerialWriter(filename); }
}

}

int main(int argc, char *argv[]) {

  int available_thread_level = -1;
  int requested_thread_level = MPI_THREAD_MULTIPLE;
  MPI_Init_thread(
      &argc, &argv, requested_thread_level, &available_thread_level);
  LOG_ERR_IF(available_thread_level != requested_thread_level,
             "Can't init MPI.");

  po::variables_map options;

  // generic options
  po::options_description generic("Generic options");

  // clang-format off
  generic.add_options()
      ("help,h", "produce this message")
      ("output,o", po::value<std::string>(), "Name of the output file.")
      ;
  // clang-format on

  // first parse cmdline and check what config file to use
  po::store(po::parse_command_line(argc, argv, generic), options);

  if (options.count("help") != 0) {
    std::cout << generic << "\n";
    std::exit(EXIT_SUCCESS);
  }

  if (options.count("output") == 0) {
    std::cout << "Missing argument `--output OUTPUT`.\n";
    std::exit(EXIT_FAILURE);
  }

  auto output_file = options["output"].as<std::string>();

  //  zisa::test_parallel_hdf5(output_file);
  //  PRINT("------------------------------------");
  zisa::test_gathered_hdf5();
  //  zisa::test_open_close_open();

  MPI_Finalize();
  return EXIT_SUCCESS;
}

#else
#include <zisa/config.hpp>
int main(int argc, char *argv[]) { LOG_ERR("Needs to be compiles with MPI"); }
#endif
