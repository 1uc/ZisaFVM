#ifndef ZISA_MPI_NUMERICAL_EXPERIMENT_HPP
#define ZISA_MPI_NUMERICAL_EXPERIMENT_HPP

#if ZISA_HAS_MPI == 1
#include <zisa/boundary/halo_exchange_bc.hpp>
#include <zisa/cli/input_parameters.hpp>
#include <zisa/io/backtrace.hpp>
#include <zisa/io/gathered_vis_info.hpp>
#include <zisa/io/gathered_visualization.hpp>
#include <zisa/io/parallel_dump_snapshot.hpp>
#include <zisa/memory/array_cell_flags.hpp>
#include <zisa/model/distributed_cfl_condition.hpp>
#include <zisa/mpi/io/hdf5_unstructured_writer.hpp>
#include <zisa/mpi/io/mpi_progress_bar.hpp>
#include <zisa/mpi/math/distributed_reference_solution.hpp>
#include <zisa/mpi/parallelization/mpi_all_reduce.hpp>
#include <zisa/mpi/parallelization/mpi_all_variables_gatherer.hpp>
#include <zisa/mpi/parallelization/mpi_halo_exchange.hpp>
#include <zisa/mpi/parallelization/mpi_single_node_array_gatherer.hpp>
#include <zisa/mpi/parallelization/mpi_single_node_array_scatterer.hpp>
#include <zisa/ode/simulation_clock.hpp>
#include <zisa/ode/time_keeper_factory.hpp>
#include <zisa/parallelization/all_variables_gatherer.hpp>
#include <zisa/parallelization/all_variables_scatterer.hpp>
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
    return this->params["io"].value("parallel_strategy", std::string("split"));
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
    auto grid = this->choose_grid();

    if (is_gathered_visualization()) {
      auto vis_info = choose_gathered_vis_info();
      auto gatherer_factory = choose_gatherer_factory();
      auto gatherer = gatherer_factory->template create_object<CellFlags, 1>();

      if (vis_info->h5_comm != MPI_COMM_NULL) {

        auto cell_flags = array<CellFlags, 1>(vis_info->n_vis_cells());
        gatherer.copy_local_patch(
            array_view(cell_flags),
            const_slice(array_const_view(grid->cell_flags),
                        0,
                        vis_info->n_local_cells));

        gatherer.receive(array_view(cell_flags));
        apply_permutation(array_view(cell_flags), *vis_info->permutation);

        std::string dirname = this->params["grid"]["file"];
        auto grid_name
            = string_format("%s/grid.h5", dirname.c_str(), mpi_comm_size);

        auto file_dims = choose_gathered_file_info();
        auto writer
            = HDF5UnstructuredWriter(grid_name, file_dims, HDF5Access::modify);

        writer.open_group("cell_flags");
        writer.unlink("ghost_cell");
        writer.close_group();

        save(writer, cell_flags, "cell_flags");
      } else {
        gatherer.send(const_slice(
            array_const_view(grid->cell_flags), 0, vis_info->n_local_cells));
      }

    } else if (is_split_visualization()) {
      auto dir = string_format("part-%04d/", mpi_rank);
      create_directory(dir);

      auto writer = HDF5SerialWriter(dir + "grid.h5");
      save(writer, *grid);

    } else if (is_unstructured_visualization()) {
      std::string dirname = this->params["grid"]["file"];
      auto grid_name
          = string_format("%s/grid.h5", dirname.c_str(), mpi_comm_size);

      auto file_dims = choose_file_dimensions();
      auto writer
          = HDF5UnstructuredWriter(grid_name, file_dims, HDF5Access::modify);

      writer.open_group("cell_flags");
      writer.unlink("ghost_cell");
      writer.close_group();
      save(writer, grid->cell_flags, "cell_flags");

    } else {
      LOG_ERR("Broken logic.");
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

    auto u1_ref = this->deduce_reference_solution(*u1);
    down_sample(u1_ref, "reference.h5");

    auto fng = this->choose_file_name_generator();
    auto steady_state_filename = fng->steady_state_filename;
    auto u_delta = load_all_vars(steady_state_filename);

    auto halo_exchange = choose_halo_exchange();
    (*halo_exchange)(*u_delta);

    for (int_t i = 0; i < u_delta->size(); ++i) {
      (*u_delta)[i] = (*u1)[i] - (*u_delta)[i];
    }

    auto grid = this->choose_grid();
    auto local_eos = this->choose_local_eos();
    auto weno_params = this->choose_weno_reference_params();
    auto local_rc_params = this->choose_local_rc_params();
    auto rc = make_reconstruction_array<NoEquilibrium,
                                        CWENO_AO,
                                        UnityScaling,
                                        typename super::eos_t,
                                        typename super::gravity_t>(
        grid, weno_params, *local_eos, this->gravity, local_rc_params);
    auto grc = std::make_shared<
        EulerGlobalReconstruction<NoEquilibrium, CWENO_AO, UnityScaling>>(
        weno_params, std::move(rc));

    auto u_delta_ref = this->deduce_reference_solution_eq(*u1, grc);
    down_sample(u_delta_ref, "delta.h5");
  }

  std::shared_ptr<AllVariables> load_all_vars(const std::string &filename) {
    if (is_gathered_visualization()) {
      auto grid = this->choose_grid();
      auto vis_info = choose_gathered_vis_info();
      auto scatter_tag = ZISA_MPI_TAG_LOAD_ALL_VARS;
      if (vis_info->h5_comm != MPI_COMM_NULL) {
        auto file_dims = choose_gathered_file_info();
        auto reader = HDF5UnstructuredReader(filename, file_dims);

        auto all_vars_patch
            = AllVariables::load(reader, all_labels<euler_var_t>());

        reverse_permutation(array_view(all_vars_patch.cvars),
                            *vis_info->permutation);

        auto darray_info
            = make_distributed_array_info(vis_info->vis_boundaries);

        auto cvars_scatterer
            = std::make_unique<MPISingleNodeArrayScatterer<double, 2>>(
                darray_info, vis_info->vis_comm, scatter_tag);

        auto avars_scatterer
            = std::make_unique<MPISingleNodeArrayScatterer<double, 2>>(
                darray_info, vis_info->vis_comm, scatter_tag + 1);

        auto scatterer = AllVariablesScatterer(std::move(cvars_scatterer),
                                               std::move(avars_scatterer),
                                               vis_info->n_local_cells,
                                               grid->n_cells);

        return std::make_shared<AllVariables>(
            scatterer.scatter(all_vars_patch));

      } else {
        auto cvars_scatterer
            = std::make_unique<MPISingleNodeArrayScatterer<double, 2>>(
                nullptr, vis_info->vis_comm, scatter_tag);

        auto avars_scatterer
            = std::make_unique<MPISingleNodeArrayScatterer<double, 2>>(
                nullptr, vis_info->vis_comm, scatter_tag + 1);

        auto scatterer = AllVariablesScatterer(std::move(cvars_scatterer),
                                               std::move(avars_scatterer),
                                               vis_info->n_local_cells,
                                               grid->n_cells);

        return std::make_shared<AllVariables>(scatterer.scatter(5, 0));
      }
    }

    if (is_unstructured_visualization()) {
      auto file_dims = choose_file_dimensions();
      auto reader = HDF5UnstructuredReader(filename, file_dims);

      return std::make_shared<AllVariables>(
          AllVariables::load(reader, all_labels<euler_var_t>()));
    }

    LOG_ERR("Implement missing case.");
  }

  void down_sample(const std::shared_ptr<ReferenceSolution> &ref_soln,
                   const std::string &filename) {

    std::vector<std::string> coarse_grid_paths
        = this->params["reference"]["coarse_grids"];

    auto serialize = [this](const std::string &filename,
                            MPI_Comm small_comm,
                            const DistributedGrid &small_dgrid,
                            const AllVariables &all_vars) {
      auto vis_info = this->compute_gathered_vis_info(small_comm, small_dgrid);
      auto gatherer_factory = this->compute_gatherer_factory(*vis_info);

      auto cvars_gatherer
          = gatherer_factory->template create_pointer<double, 2>();

      // FIXME missing adv variables.

      auto n_local_cells = vis_info->n_local_cells;
      auto n_vis_cells = vis_info->n_vis_cells();
      auto all_var_gatherer = AllVariablesGatherer(
          std::move(cvars_gatherer), nullptr, n_local_cells, n_vis_cells);

      auto gathered_all_vars = all_var_gatherer.gather(all_vars);

      if (vis_info->h5_comm != MPI_COMM_NULL) {
        auto array_info = compute_gathered_file_info(*vis_info);
        auto writer = HDF5UnstructuredWriter(filename, array_info);

        double t_end = 0.0;
        int_t n_steps = 0;

        apply_permutation(array_view(gathered_all_vars.cvars),
                          *vis_info->permutation);

        auto labels = all_labels<euler_var_t>();
        save_state(writer, gathered_all_vars, t_end, n_steps, labels);
      }
    };

    auto fine_grid = this->choose_grid();
    auto mask = super::boundary_mask();

    auto interpolation = [&ref_soln](int_t i, const XYZ &x, int_t k) {
      return ref_soln->q_ref(x, k, i);
    };

    auto dref = DistributedReferenceSolution(
        serialize, fine_grid, interpolation, euler_var_t::size());

    for (const auto &coarse_base_folder : coarse_grid_paths) {
      auto coarse_folder
          = std::filesystem::path(coarse_base_folder) / "partitioned";
      dref.compute_and_save(
          filename, coarse_folder, small_comm_size(coarse_folder), mask);
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
    return string_format("%s/partitioned/%d/subgrid-%04d.msh.h5",
                         dirname.c_str(),
                         mpi_comm_size,
                         mpi_rank);
  }

  std::shared_ptr<Grid> compute_grid() const override {
    // 1. Load in a oversized chunk.
    // 2. Decide how much we really need based on the stencil.
    // 3. Generate correctly sized grid.

    // Side effect: store the stencils.

    auto quad_deg = this->choose_volume_deg();
    auto super_subgrid = zisa::load_grid(subgrid_name(), quad_deg);
    auto super_sub_dgrid = zisa::load_distributed_grid(subgrid_name());

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
    if (is_gathered_visualization()) {
      return compute_gathered_visualization();
    }

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
    auto local_eos = this->compute_local_eos();

    // TODO here we just made this only work for Euler.
    return std::make_shared<ParallelDumpSnapshot<typename super::eos_t>>(
        local_eos, fng, file_dims);
  }

  std::shared_ptr<GatheredVisInfo> choose_gathered_vis_info() {
    if (gathered_vis_info_ == nullptr) {
      gathered_vis_info_ = compute_gathered_vis_info();
    }

    return gathered_vis_info_;
  }

  std::shared_ptr<GatheredVisInfo> compute_gathered_vis_info() {
    auto dgrid = choose_distributed_grid();
    return compute_gathered_vis_info(mpi_comm, *dgrid);
  }

  std::shared_ptr<GatheredVisInfo>
  compute_gathered_vis_info(MPI_Comm world_comm, const DistributedGrid &dgrid) {
    auto world_rank = mpi_rank;
    auto world_size = mpi_comm_size;

    auto n_local_cells = int_t(std::count(
        dgrid.partition.begin(), dgrid.partition.end(), int_t(world_rank)));

    auto n_vis_pes
        = zisa::min(int(this->params["io"]["n_writers"]), world_size);
    auto pes_per_block = (world_size + n_vis_pes - 1) / n_vis_pes;
    auto vis_id = world_rank / pes_per_block;

    auto vis_comm = zisa::mpi::comm_split(world_comm, vis_id, world_rank);
    auto vis_size = zisa::mpi::size(vis_comm);
    auto vis_rank = zisa::mpi::rank(vis_comm);
    auto vis_tag = ZISA_MPI_TAG_GATHERED_VIS;

    auto h5_comm = zisa::mpi::comm_split(
        world_comm, vis_rank == 0 ? 0 : MPI_UNDEFINED, world_rank);

    {
      std::string name = "VisComm-";
      name += std::to_string(vis_id);
      MPI_Comm_set_name(vis_comm, name.c_str());
    }

    if (h5_comm != MPI_COMM_NULL) {
      auto vis_cells_per_pe = array<int_t, 1>(shape_t<1>{vis_size});
      vis_cells_per_pe[0] = n_local_cells;
      zisa::mpi::gather(array_view(vis_cells_per_pe), 0, vis_comm);
      auto n_vis_cells = std::accumulate(
          vis_cells_per_pe.begin(), vis_cells_per_pe.end(), int_t(0));

      auto gids = array<int_t, 1>(shape_t<1>{n_vis_cells});
      auto vis_boundaries = array<int_t, 1>(shape_t<1>{vis_size + 1});
      vis_boundaries[0] = 0;
      for (int i = 0; i < vis_size; ++i) {
        vis_boundaries[i + 1] = vis_boundaries[i] + vis_cells_per_pe[i];
      }

      for (int_t i = 0; i < n_local_cells; ++i) {
        gids[i] = dgrid.global_cell_indices[i];
      }

      for (int p = 1; p < vis_size; ++p) {
        auto i0 = vis_boundaries[p];
        auto i1 = vis_boundaries[p + 1];
        zisa::mpi::recv(slice(array_view(gids), i0, i1), p, vis_tag, vis_comm);
      }

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

      return std::make_shared<GatheredVisInfo>(
          GatheredVisInfo{std::move(vis_boundaries),
                          std::move(gids),
                          sigma,
                          n_local_cells,
                          vis_comm,
                          h5_comm});

    } else {
      zisa::mpi::gather(array_view(shape_t<1>{1}, &n_local_cells), 0, vis_comm);
      zisa::mpi::send(const_slice(array_const_view(dgrid.global_cell_indices),
                                  0,
                                  n_local_cells),
                      0,
                      vis_tag,
                      vis_comm);

      return std::make_shared<GatheredVisInfo>(
          GatheredVisInfo{array<int_t, 1>(),
                          array<int_t, 1>(),
                          nullptr,
                          n_local_cells,
                          vis_comm,
                          h5_comm});
    }
  }

  std::shared_ptr<HDF5UnstructuredFileDimensions> choose_gathered_file_info() {
    if (gathered_file_info_ == nullptr) {
      gathered_file_info_ = compute_gathered_file_info();
    }

    return gathered_file_info_;
  }

  std::shared_ptr<HDF5UnstructuredFileDimensions> compute_gathered_file_info() {
    auto vis_info = choose_gathered_vis_info();
    return compute_gathered_file_info(*vis_info);
  }

  std::shared_ptr<HDF5UnstructuredFileDimensions>
  compute_gathered_file_info(const GatheredVisInfo &vis_info) {

    if (vis_info.h5_comm != MPI_COMM_NULL) {
      const auto &gids = vis_info.vis_file_ids;
      std::vector<hsize_t> hids(gids.size());
      for (int_t i = 0; i < hids.size(); ++i) {
        hids[i] = integer_cast<hsize_t>(gids[i]);
      }

      return make_hdf5_unstructured_file_dimensions(
          vis_info.n_vis_cells(), hids, vis_info.h5_comm);
    } else {
      return nullptr;
    }
  }

  std::shared_ptr<MPISingleNodeArrayGathererFactory> choose_gatherer_factory() {
    if (gatherer_factory_ == nullptr) {
      gatherer_factory_ = compute_gatherer_factory();
    }
    return gatherer_factory_;
  }

  std::shared_ptr<MPISingleNodeArrayGathererFactory>
  compute_gatherer_factory() {
    auto vis_info = choose_gathered_vis_info();
    return compute_gatherer_factory(*vis_info);
  }

  std::shared_ptr<MPISingleNodeArrayGathererFactory>
  compute_gatherer_factory(const GatheredVisInfo &vis_info) {
    auto vis_comm = vis_info.vis_comm;
    auto vis_tag = ZISA_MPI_TAG_GATHERED_VIS;

    auto h5_comm = vis_info.h5_comm;

    if (h5_comm != MPI_COMM_NULL) {
      auto darray_info = make_distributed_array_info(vis_info.vis_boundaries);
      return std::make_shared<MPISingleNodeArrayGathererFactory>(
          darray_info, vis_comm, vis_tag + 1);
    } else {
      return std::make_shared<MPISingleNodeArrayGathererFactory>(
          nullptr, vis_comm, vis_tag + 1);
    }
  }

  std::shared_ptr<Visualization> compute_gathered_visualization() {
    auto vis_info = choose_gathered_vis_info();
    auto n_local_cells = vis_info->n_local_cells;
    auto n_vis_cells = vis_info->n_vis_cells();
    auto h5_comm = vis_info->h5_comm;

    auto gatherer_factory = choose_gatherer_factory();

    auto cvars_gatherer
        = gatherer_factory->template create_pointer<double, 2>();

    auto avars_gatherer
        = gatherer_factory->template create_pointer<double, 2>();
    auto all_var_gatherer
        = std::make_unique<AllVariablesGatherer>(std::move(cvars_gatherer),
                                                 std::move(avars_gatherer),
                                                 n_local_cells,
                                                 n_vis_cells);

    if (h5_comm != MPI_COMM_NULL) {
      auto fng = this->choose_file_name_generator();

      auto file_dims = choose_gathered_file_info();
      auto local_eos = this->compute_local_eos(file_dims->n_cells_local);
      auto dump_snapshot
          = std::make_shared<ParallelDumpSnapshot<typename super::eos_t>>(
              local_eos, fng, file_dims);

      auto all_var_dims = this->choose_all_variable_dims();
      all_var_dims.n_cells = n_vis_cells;

      return std::make_shared<GatheredVisualization>(
          std::move(all_var_gatherer),
          vis_info->permutation,
          dump_snapshot,
          all_var_dims);

    } else {
      return std::make_shared<GatheredVisualization>(
          std::move(all_var_gatherer),
          nullptr,
          nullptr,
          AllVariablesDimensions{0, 0, 0});
    }
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
  mutable std::shared_ptr<GatheredVisInfo> gathered_vis_info_ = nullptr;
  mutable std::shared_ptr<HDF5UnstructuredFileDimensions> gathered_file_info_
      = nullptr;
  mutable std::shared_ptr<MPISingleNodeArrayGathererFactory> gatherer_factory_
      = nullptr;
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
