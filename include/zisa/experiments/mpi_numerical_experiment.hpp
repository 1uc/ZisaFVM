#ifndef ZISA_MPI_NUMERICAL_EXPERIMENT_HPP
#define ZISA_MPI_NUMERICAL_EXPERIMENT_HPP

#if ZISA_HAS_MPI == 1
#include <zisa/boundary/halo_exchange_bc.hpp>
#include <zisa/cli/input_parameters.hpp>
#include <zisa/io/gathered_visualization.hpp>
#include <zisa/io/parallel_dump_snapshot.hpp>
#include <zisa/model/distributed_cfl_condition.hpp>
#include <zisa/mpi/io/hdf5_unstructured_writer.hpp>
#include <zisa/mpi/io/mpi_progress_bar.hpp>
#include <zisa/mpi/math/distributed_reference_solution.hpp>
#include <zisa/mpi/parallelization/mpi_all_reduce.hpp>
#include <zisa/mpi/parallelization/mpi_all_variables_gatherer.hpp>
#include <zisa/mpi/parallelization/mpi_halo_exchange.hpp>
#include <zisa/mpi/parallelization/mpi_single_node_array_gatherer.hpp>
#include <zisa/ode/simulation_clock.hpp>
#include <zisa/ode/time_keeper_factory.hpp>
#include <zisa/parallelization/all_variables_gatherer.hpp>
#include <zisa/parallelization/distributed_grid.hpp>
#include <zisa/parallelization/domain_decomposition.hpp>

namespace zisa {

template <class NMExp>
class MPINumericalExperiment : public NMExp {
private:
  using super = NMExp;

public:
  explicit MPINumericalExperiment(const InputParameters &params)
      : super(params),
        mpi_rank(zisa::mpi::rank(mpi_comm)),
        mpi_comm_size(zisa::mpi::size(mpi_comm)) {}

protected:
  GMSHElementType deduce_element_type(int_t max_neighbours) const {
    return max_neighbours == 3 ? GMSHElementType::triangle
                               : GMSHElementType::tetrahedron;
  }

  std::string parallel_visualization_strategy() const {
    return this->params["io"].value("parallel_strategy",
                                    std::string("unstructured"));
  }

  bool is_unstructured_visualization() const {
    return parallel_visualization_strategy() == "unstructured";
  }

  bool is_gathered_visualization() const {
    return parallel_visualization_strategy() == "gathered";
  }

  bool is_split_visualization() const {
    return parallel_visualization_strategy() == "split";
  }

  void write_grid() override {
    if (is_gathered_visualization()) {
      LOG_WARN(
          "This needs to be reimplmented. In the meantime, simply copy over.");
    } else if (is_split_visualization()) {
      auto dir = string_format("part-%04d/", mpi_rank);
      create_directory(dir);

      auto writer = HDF5SerialWriter(dir + "grid.h5");
      save(writer, *this->grid_);
    }
  }

  void print_grid_info() override {
    if (mpi_rank == 0) {
      super::print_grid_info();
    }
  }

  void do_post_run(const std::shared_ptr<AllVariables> &u1) override {
    if (!has_key(this->params, "reference")) {
      // Post processing the reference solution is not requested.
      return;
    }

    auto u1_ref = this->deduce_reference_solution(u1);
    down_sample(u1_ref, "reference.h5");

    auto fng = this->choose_file_name_generator();
    auto file_dims = choose_file_dimensions();
    auto steady_state_filename = fng->steady_state_filename;
    auto reader = HDF5UnstructuredReader(steady_state_filename, file_dims);
    auto u_delta = std::make_shared<AllVariables>(
        AllVariables::load(reader, all_labels<euler_var_t>()));

    for (int_t i = 0; i < u_delta->size(); ++i) {
      (*u_delta)[i] = (*u1)[i] - (*u_delta)[i];
    }

    auto u_delta_ref = this->deduce_reference_solution_eq(
        u_delta, NoEquilibrium{}, UnityScaling{});
    down_sample(u_delta_ref, "delta.h5");
  }

  void down_sample(const std::shared_ptr<ReferenceSolution> &ref_soln,
                   const std::string &filename) {

    std::vector<std::string> coarse_grid_paths
        = this->params["reference"]["coarse_grids"];

    auto serialize = [this](HDF5Writer &writer, const AllVariables &all_vars) {
      double t_end = 0.0;
      int_t n_steps = 0;
      save_full_state(writer, *this->euler, all_vars, t_end, n_steps);
    };

    auto fine_grid = this->choose_grid();
    auto mask = this->boundary_mask();

    auto interpolation = [&ref_soln](int_t i, const XYZ &x, int_t k) {
      return ref_soln->q_ref(x, k, i);
    };

    auto dref = DistributedReferenceSolution(
        serialize, fine_grid, interpolation, euler_var_t::size());

    for (const auto &coarse_grid_folder : coarse_grid_paths) {
      dref.compute_and_save(filename,
                            coarse_grid_folder,
                            small_comm_size(coarse_grid_folder),
                            mask);
    }
  }

  int small_comm_size(const std::string &coarse_folder) {
    auto sizes = std::vector<int>();
    for (auto &p : std::filesystem::directory_iterator(coarse_folder)) {
      if (p.is_directory()) {
        auto b = p.path().filename().string();
        auto s = stoi(b);
        if (string_format("%d", s) == b) {
          if (s <= mpi_comm_size) {
            sizes.push_back(s);
          }
        }
      }
    }

    LOG_ERR_IF(sizes.empty(), "Failed to deduce any valid `small_comm_size`.");
    return *std::max_element(sizes.begin(), sizes.end());
  }

  void do_post_process() override {}

  std::shared_ptr<array<StencilFamily, 1>> choose_stencils() const override {
    LOG_ERR_IF(this->stencils_ == nullptr,
               "Trying to use stencils before they were computed.");
    return this->stencils_;
  }

  std::shared_ptr<DistributedGrid> choose_distributed_grid() const {
    LOG_ERR_IF(this->distributed_grid_ == nullptr,
               "Trying to use stencils before they were computed.");
    return this->distributed_grid_;
  }

  std::shared_ptr<array<StencilFamily, 1>>
  compute_stencils(const Grid &grid) const override {
    auto stencil_params = this->choose_stencil_params();

    return std::make_shared<array<StencilFamily, 1>>(
        compute_stencil_families(grid, stencil_params));
  }

  std::string subgrid_name() const {
    std::string dirname = this->params["grid"]["file"];
    return string_format(
        "%s/%d/subgrid-%04d.msh.h5", dirname.c_str(), mpi_comm_size, mpi_rank);
  }

  std::shared_ptr<Grid> compute_grid() const override {
    // 1. Load in a oversized chunk.
    // 2. Decide how much we really need based on the stencil.
    // 3. Generate correctly sized grid.

    // Side effect: store the stencils.

    auto super_subgrid = zisa::load_grid(subgrid_name());
    auto super_sub_dgrid = zisa::load_distributed_grid(subgrid_name());

    // FIXME this need to be extracted & stored.
    auto stencils = compute_stencils(*super_subgrid);

    auto is_interior = [this, &partition = super_sub_dgrid.partition](int_t i) {
      return partition[i] == integer_cast<int_t>(mpi_rank);
    };

    auto is_needed
        = StencilBasedIndicator(*super_subgrid, *stencils, is_interior);

    auto [local_vertex_indices, local_vertices, super_sub_indices]
        = extract_subgrid_v2(*super_subgrid, is_needed);

    this->stencils_ = std::make_shared<array<StencilFamily, 1>>(
        extract_stencils(*stencils,
                         super_subgrid->neighbours,
                         is_interior,
                         super_sub_indices));

    distributed_grid_ = std::make_shared<DistributedGrid>(
        extract_distributed_subgrid(super_sub_dgrid, super_sub_indices));

    return compute_local_grid(std::move(local_vertex_indices),
                              std::move(local_vertices));
  }

  std::shared_ptr<Grid> compute_local_grid(array<int_t, 2> vertex_indices,
                                           array<XYZ, 1> vertices) const {

    auto quad_deg = this->choose_volume_deg();
    auto max_neighbours = vertex_indices.shape(1);

    auto grid = std::make_shared<Grid>(deduce_element_type(max_neighbours),
                                       std::move(vertices),
                                       std::move(vertex_indices),
                                       quad_deg);

    this->enforce_cell_flags(*grid);
    return grid;
  }

  std::function<bool(const Grid &, int_t i)> boundary_mask() const override {
    auto physical_mask = super::boundary_mask();
    auto dgrid = choose_distributed_grid();
    return [physical_mask, dgrid, mpi_rank = this->mpi_rank](const Grid &grid,
                                                             int_t i) {
      return dgrid->partition[i] != integer_cast<int_t>(mpi_rank)
             || physical_mask(grid, i);
    };
  }

  std::shared_ptr<AllVariables> load_initial_conditions() override {
    LOG_ERR("Implement first.");
  }

  std::shared_ptr<StepRejection> choose_step_rejection() override {
    // Only RejectNothing works without implementing more.
    return std::dynamic_pointer_cast<RejectNothing>(
        super::choose_step_rejection());
  }

  std::shared_ptr<Visualization> compute_visualization() override {
    if (is_unstructured_visualization()) {
      return compute_unstructured_visualization();
    }

    if (is_split_visualization()) {
      return super::compute_visualization();
    }

    LOG_ERR("Invalid config for visualization.");
  }

  std::shared_ptr<HDF5UnstructuredFileDimensions>
  choose_file_dimensions() const {
    return compute_file_dimensions();
  }

  std::shared_ptr<HDF5UnstructuredFileDimensions>
  compute_file_dimensions() const {
    auto dgrid = choose_distributed_grid();
    auto global_ids = std::vector<hsize_t>();

    for (int_t i = 0; i < dgrid->global_cell_indices.size(); ++i) {
      if (dgrid->partition[i] == integer_cast<int_t>(mpi_rank)) {
        global_ids.push_back(dgrid->global_cell_indices[i]);
      }
    }

    return make_hdf5_unstructured_file_dimensions(
        dgrid->partition.shape(0), global_ids, mpi_comm);
  }

  std::shared_ptr<Visualization> compute_unstructured_visualization() {
    const auto &fng = this->choose_file_name_generator();
    auto file_dims = choose_file_dimensions();

    // TODO here we just made this only work for Euler.
    return std::make_shared<ParallelDumpSnapshot<typename super::euler_t>>(
        this->euler, fng, file_dims);
  }

  std::shared_ptr<SimulationClock> choose_simulation_clock() override {
    auto time_keeper_params = TimeKeeperParameters(this->params);
    auto plotting_params = PlottingStepsParameters(this->params);

    auto time_keeper = make_time_keeper(time_keeper_params);
    auto plotting_steps
        = make_plotting_steps(plotting_params, time_keeper_params);

    return std::make_shared<MPISimulationClock>(
        time_keeper, plotting_steps, mpi_comm);
  }

  std::shared_ptr<CFLCondition> choose_cfl_condition() override {
    auto cfl = super::choose_cfl_condition();

    auto op = ReductionOperation::min;
    auto all_reduce = std::make_shared<MPIAllReduce>(op, mpi_comm);

    return std::make_shared<DistributedCFLCondition>(cfl, all_reduce);
  }

  std::shared_ptr<HaloExchange> choose_halo_exchange() {
    if (halo_exchange_ == nullptr) {
      halo_exchange_ = compute_halo_exchange();
    }

    return halo_exchange_;
  }

  std::shared_ptr<HaloExchange> compute_halo_exchange() {
    auto i_need_this = std::vector<HaloRemoteInfo>();
    auto local_halo_info = std::vector<HaloReceiveInfo>();

    auto dgrid = choose_distributed_grid();
    const auto &partition = dgrid->partition;
    const auto &global_indices = dgrid->global_cell_indices;

    auto n_owned_cells
        = std::count(partition.begin(), partition.end(), mpi_rank);
    auto n_cells = partition.shape(0);

    int_t i_start = integer_cast<int_t>(n_owned_cells);
    while (i_start < n_cells) {
      auto p = integer_cast<int>(partition(i_start));

      int_t i_end = i_start;
      while (i_end < n_cells && integer_cast<int>(partition[i_end]) == p) {
        ++i_end;
      }

      auto n_patch = i_end - i_start;
      auto indices = array<int_t, 1>(n_patch);
      for (int_t i = i_start; i < i_end; ++i) {
        indices[i - i_start] = global_indices[i];
      }

      i_need_this.emplace_back(p, indices);
      local_halo_info.emplace_back(p, i_start, i_end);

      i_start = i_end;
    }

    return std::make_shared<MPIHaloExchange>(make_mpi_halo_exchange(
        *dgrid,
        {std::move(i_need_this), std::move(local_halo_info)},
        mpi_comm));
  }

  std::shared_ptr<BoundaryCondition> compute_boundary_condition() override {
    auto bc = super::compute_boundary_condition();
    auto halo_exchange = choose_halo_exchange();
    return std::make_shared<HaloExchangeBC>(bc, halo_exchange);
  }

  std::shared_ptr<ProgressBar> choose_progress_bar() override {
    auto serial_bar = super::choose_progress_bar();
    return std::make_shared<MPIProgressBar>(serial_bar, mpi_comm);
  }

  std::shared_ptr<FileNameGenerator> compute_file_name_generator() override {
    const auto &fn_params = this->params["io"]["filename"];

    std::string stem = fn_params["stem"];
    std::string dir
        = (is_split_visualization() ? string_format("part-%04d/", mpi_rank)
                                    : std::string("./"));

    return make_file_name_generator(
        dir, stem, fn_params["pattern"], fn_params["suffix"]);
  }

  void write_debug_output() override {
    if (mpi_rank == 0) {
      super::write_debug_output();
    }
  }

protected:
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  int mpi_rank;
  int mpi_comm_size;

  mutable std::shared_ptr<DistributedGrid> distributed_grid_ = nullptr;
  mutable std::shared_ptr<HaloExchange> halo_exchange_ = nullptr;
};
}
#else
#include <zisa/cli/input_parameters.hpp>

namespace zisa {

template <class NMExp>
class MPINumericalExperiment : public NMExp {
private:
  using super = NMExp;

public:
  explicit MPINumericalExperiment(const InputParameters &) {
    LOG_ERR("MPI requested, but not build with `-DZISA_HAS_MPI=1`.");
  }
};
}

#endif
#endif // ZISA_MPI_NUMERICAL_EXPERIMENT_HPP
