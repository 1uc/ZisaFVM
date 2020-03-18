#ifndef ZISA_MPI_NUMERICAL_EXPERIMENT_HPP
#define ZISA_MPI_NUMERICAL_EXPERIMENT_HPP

#if ZISA_HAS_MPI == 1
#include <zisa/boundary/halo_exchange_bc.hpp>
#include <zisa/cli/input_parameters.hpp>
#include <zisa/io/gathered_visualization.hpp>
#include <zisa/io/mpi_progress_bar.hpp>
#include <zisa/model/distributed_cfl_condition.hpp>
#include <zisa/ode/simulation_clock.hpp>
#include <zisa/ode/time_keeper_factory.hpp>
#include <zisa/parallelization/all_variables_gatherer.hpp>
#include <zisa/parallelization/domain_decomposition.hpp>
#include <zisa/parallelization/mpi_all_reduce.hpp>
#include <zisa/parallelization/mpi_all_variables_gatherer.hpp>
#include <zisa/parallelization/mpi_halo_exchange.hpp>
#include <zisa/parallelization/mpi_single_node_array_gatherer.hpp>

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
                                    std::string("gathered"));
  }

  bool is_gathered_visualization() const {
    return parallel_visualization_strategy() == "gathered";
  }

  bool is_split_visualization() const {
    return parallel_visualization_strategy() == "split";
  }

  void write_grid() override {
    if (is_gathered_visualization()) {
      if (mpi_rank == 0) {
        auto writer = HDF5SerialWriter("grid.h5");
        save(writer, *this->full_grid_);
      }
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
    // Check to see if it needs gathering.
    if (u1->dims().n_cells == this->full_grid_->n_cells) {
      super::do_post_process();
    } else {
      const auto &boundaries = partitioned_grid_->boundaries;
      auto gatherer = make_mpi_all_variables_gatherer(
          ZISA_MPI_TAG_ALL_VARS_POSTPROCESS, mpi_comm, boundaries);

      auto full_all_vars
          = std::make_shared<AllVariables>(gatherer->gather(*u1));

      if (mpi_rank == 0) {
        // Restore original order.
        const auto &sigma = factor_permutation(partitioned_grid_->permutation);
        reverse_permutation(array_view(full_all_vars->cvars), sigma);
        reverse_permutation(array_view(full_all_vars->avars), sigma);

        super::do_post_run(full_all_vars);
      }
    }
  }

  void do_post_process() override {
    if (mpi_rank == 0) {
      super::do_post_process();
    }
  }

  std::shared_ptr<Grid> compute_full_grid() const override {
    int_t quad_deg = this->choose_volume_deg();

    if (mpi_rank == 0) {
      auto vertex_indices_gbl = array<int_t, 2>{};
      auto vertices_gbl = array<XYZ, 1>{};
      {
        auto reader = HDF5SerialReader(this->params["grid"]["file"]);
        vertex_indices_gbl = array<int_t, 2>::load(reader, "vertex_indices");
        vertices_gbl = array<XYZ, 1>::load(reader, "vertices");
      }

      auto n_cells = vertex_indices_gbl.shape(0);
      auto max_neighbours = vertex_indices_gbl.shape(1);
      auto n_vertices = vertices_gbl.shape(0);

      auto sizes = shape_t<3>{n_cells, n_vertices, max_neighbours};
      zisa::mpi::bcast(array_view(sizes.shape(), sizes.raw()), 0, mpi_comm);
      zisa::mpi::bcast(array_view(vertex_indices_gbl), 0, mpi_comm);
      zisa::mpi::bcast(array_view(vertices_gbl), 0, mpi_comm);

      auto grid = std::make_shared<Grid>(deduce_element_type(max_neighbours),
                                    std::move(vertices_gbl),
                                    std::move(vertex_indices_gbl),
                                    quad_deg);
      this->enforce_cell_flags(*grid);
      return grid;

    } else {
      auto sizes = shape_t<3>{};
      mpi::bcast(array_view(sizes.shape(), sizes.raw()), 0, mpi_comm);
      auto n_cells = sizes[0];
      auto n_vertices = sizes[1];
      auto max_neighbours = sizes[2];

      auto vertex_indices_gbl = array<int_t, 2>({n_cells, max_neighbours});
      auto vertices_gbl = array<XYZ, 1>(n_vertices);

      zisa::mpi::bcast(array_view(vertex_indices_gbl), 0, mpi_comm);
      zisa::mpi::bcast(array_view(vertices_gbl), 0, mpi_comm);

      auto grid = std::make_shared<Grid>(deduce_element_type(max_neighbours),
                                    std::move(vertices_gbl),
                                    std::move(vertex_indices_gbl),
                                    quad_deg);
      this->enforce_cell_flags(*grid);
      return grid;
    }
  }

  std::shared_ptr<PartitionedGrid>
  partition_grid(const Grid &full_grid,
                 const array<StencilFamily, 1> &stencils) const {

    if (mpi_rank == 0) {
      partitioned_grid_ = std::make_shared<PartitionedGrid>(
          compute_partitioned_grid(full_grid, stencils, int_t(mpi_comm_size)));
    } else {
      auto n_cells = full_grid.n_cells;
      partitioned_grid_ = std::make_shared<PartitionedGrid>(
          array<int_t, 1>(n_cells),
          array<int_t, 1>(mpi_comm_size + 1),
          array<int_t, 1>(n_cells));
    }

    zisa::mpi::bcast(array_view(partitioned_grid_->partition), 0, mpi_comm);
    zisa::mpi::bcast(array_view(partitioned_grid_->boundaries), 0, mpi_comm);
    zisa::mpi::bcast(array_view(partitioned_grid_->permutation), 0, mpi_comm);

    return partitioned_grid_;
  }

  std::shared_ptr<array<StencilFamily, 1>> choose_stencils() const override {
    LOG_ERR_IF(this->stencils_ == nullptr,
               "Trying to use stencils before they were computed.");
    return this->stencils_;
  }

  std::shared_ptr<array<StencilFamily, 1>>
  compute_stencils(const Grid &grid) const override {
    auto stencil_params = this->choose_stencil_params();

    return std::make_shared<array<StencilFamily, 1>>(
        compute_stencil_families(grid, stencil_params));
  }

  std::shared_ptr<Grid> compute_full_renumbered_grid(const Grid &grid) const {
    const auto &vertices = grid.vertices;
    const auto &vertex_indices = grid.vertex_indices;
    auto element_type = deduce_element_type(grid.max_neighbours);
    auto quad_deg = this->choose_volume_deg();

    return std::make_shared<Grid>(
        element_type,
        vertices,
        renumbered_vertex_indices(vertex_indices,
                                  partitioned_grid_->permutation),
        quad_deg);
  }

  std::shared_ptr<Grid> compute_grid() const override {
    auto full_grid = this->choose_full_grid();
    auto stencils = compute_stencils(*full_grid);

    auto partitioned_grid = partition_grid(*full_grid, *stencils);

    auto [local_vertex_indices, local_vertices, local_stencils, halo]
        = extract_subgrid(*full_grid,
                          *partitioned_grid,
                          *stencils,
                          integer_cast<int_t>(mpi_rank));

    halo_exchange_ = std::make_shared<MPIHaloExchange>(
        make_mpi_halo_exchange(halo, mpi_comm));

    this->full_stencils_ = stencils;
    this->stencils_
        = std::make_shared<array<StencilFamily, 1>>(std::move(local_stencils));

    return local_grid(std::move(local_vertex_indices),
                      std::move(local_vertices));
  }

  std::shared_ptr<Grid> local_grid(array<int_t, 2> vertex_indices,
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

    if (is_split_visualization()) {
      return super::compute_visualization();
    }

    LOG_ERR("Invalid config for visualization.");
  }

  std::shared_ptr<Visualization> compute_gathered_visualization() {
    auto single_node_vis = super::compute_visualization();

    const auto &boundaries = partitioned_grid_->boundaries;
    auto gatherer = make_mpi_all_variables_gatherer(
        ZISA_MPI_TAG_ALL_VARS_VISUALIZATION, mpi_comm, boundaries);

    auto all_var_dims = this->choose_all_variable_dims();
    all_var_dims.n_cells = partitioned_grid_->boundaries[mpi_comm_size];

    auto permutation = std::make_shared<Permutation>(
        factor_permutation(partitioned_grid_->permutation));

    return std::make_shared<GatheredVisualization>(std::move(gatherer),
                                                   std::move(permutation),
                                                   single_node_vis,
                                                   all_var_dims);
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

  std::shared_ptr<BoundaryCondition> compute_boundary_condition() override {
    auto bc = super::compute_boundary_condition();
    return std::make_shared<HaloExchangeBC>(bc, halo_exchange_);
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

  mutable std::shared_ptr<PartitionedGrid> partitioned_grid_ = nullptr;
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
