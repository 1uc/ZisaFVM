#ifndef ZISA_MPI_NUMERICAL_EXPERIMENT_HPP
#define ZISA_MPI_NUMERICAL_EXPERIMENT_HPP

#include <zisa/cli/input_parameters.hpp>
#include <zisa/model/distributed_cfl_condition.hpp>
#include <zisa/parallelization/all_reduce.hpp>
#include <zisa/parallelization/all_variables_gatherer.hpp>
#include <zisa/parallelization/mpi_halo_exchange.hpp>
#include <zisa/parallelization/mpi_single_node_array_gatherer.hpp>

namespace zisa {

template <class NMExp>
class MPINumericalExperiment : public NMExp {
private:
  using super = NMExp;

public:
  MPINumericalExperiment(const InputParameters &params)
      : super(params),
        mpi_comm(MPI_COMM_WORLD),
        mpi_rank(zisa::mpi::rank(mpi_comm)),
        mpi_comm_size(zisa::mpi::size(mpi_comm)) {}

protected:
  void write_grid() const override {
    if (mpi_rank == 0) {
      auto writer = HDF5SerialWriter("grid.h5");
      save(writer, *full_grid_);
    }
  }

  void print_grid_info() override {
    if (mpi_rank == 0) {
      super::print_grid_info();
    }
  }

  void do_post_run(const std::shared_ptr<AllVariables> &u1) {
    LOG_ERR("Implement first.");
  }

  void do_post_process(const std::shared_ptr<AllVariables> &u1) {
    LOG_ERR("Implement first.");
  }

  Grid load_full_grid() {
    int_t quad_deg = choose_volume_deg();

    if (rank == 0) {
      auto gmsh_data = GMSHData(params["grid"]["file"]);
      auto &vertex_indices_gbl = gmsh_data.vertex_indices;
      auto &vertices_gbl = gmsh_data.vertices;

      auto n_cells = vertex_indices_gbl.shape(0);
      auto max_neighbours = vertex_indices_gbl.shape(1);
      auto n_vertices = vertices.shape(0);

      auto sizes = shape_t<3>{n_cells, n_vertices, max_neighbours};
      mpi::bcast(sizes, 0, mpi_comm);

      zisa::mpi::bcast(vertex_indices_gbl, 0, mpi_comm);
      zisa::mpi::bcast(vertices_gbl, 0, mpi_comm);

      return Grid(
          std::move(vertex_indices_gbl), std::move(vertices_gbl), quad_deg);
    } else {
      auto sizes = shape_t<3>{n_cells, n_vertices, max_neighbours};
      mpi::bcast(sizes, 0, mpi_comm);

      auto [n_cells, n_vertices, max_neighbours] = sizes;

      auto vertices_gbl = array<XYZ, 1>(n_vertices);
      auto vertex_indices_gbl = array<int_t, 2>({n_cells, max_neighbours});

      zisa::mpi::bcast(vertex_indices_gbl, 0, mpi_comm);
      zisa::mpi::bcast(vertices_gbl, 0, mpi_comm);

      return Grid(
          std::move(vertex_indices_gbl), std::move(vertices_gbl), quad_deg);
    }
  };

  std::shared_ptr<PartitionedGrid>
  partition_grid(const Grid &full_grid,
                 const std::vector<StencilFamily> &stencils) {

    if (mpi_rank == 0) {
      partitioned_grid_ = std::make_shared<PartitionedGrid>(
          compute_partitioned_grid(full_grid, stencils, n_parts));
    } else {
      auto n_cells = full_grid.n_cells;
      partitioned_grid_
          = std::make_shared<PartitionedGrid>(array<int_t, 1>(n_cells),
                                              array<int_t, 1>(mpi_comm_size),
                                              array<int_t, 1>(n_cells));
    }

    zisa::mpi::bcast(partitioned_grid_->partition, 0, mpi_comm);
    zisa::mpi::bcast(partitioned_grid_->boundaries, 0, mpi_comm);
    zisa::mpi::bcast(partitioned_grid_->permutation, 0, mpi_comm);

    return partitioned_grid;
  }

  void compute_grid() override {
    auto full_grid = load_full_grid();
    auto stencils = compute_stencils(full_grid);
    auto partitioned_grid = partition_grid(full_grid);

    auto [local_vertex_indices, local_vertices, local_stencils, halo]
        = extract_subgrid(full_grid, partitioned_grid, stencils, mpi_rank);

    halo_exchange_ = make_mpi_halo_exchange(halo, mpi_comm);
    stencils_ = local_stencils;

    return local_grid(std::move(local_vertex_indices),
                      std::move(local_vertices));
  }

  std::shared_ptr<Grid> local_grid(array<int_t, 2> local_vertices,
                                   array<XYZ, 1> vertices) {

    auto quad_deg = choose_volume_deg();
    auto element_type = grid.max_neighbours == 3
                            ? zisa::GMSHElementType::triangle
                            : zisa::GMSHElementType::tetrahedron;

    return std::make_shared<Grid>(element_type,
                                  std::move(local_vertices),
                                  std::move(vertices),
                                  element_type,
                                  quad_deg);
  }

  std::shared_ptr<AllVariables> load_initial_conditions() override {
    LOG_ERR("Implement first.");
  }

  std::shared_ptr<Visualization> choose_visualization() override {
    auto single_node_vis = super::choose_visualization();

    auto array_info
        = std::make_shared<DistributedArrayInfo>(partitioned_grid_->boundaries);

    auto cvars_gatherer
        = std::make_unique<MPISingleNodeArrayGatherer<double, 2>>(
            array_info, mpi_comm, 388932);

    auto avars_gatherer
        = std::make_unique<MPISingleNodeArrayGatherer<double, 2>>(
            array_info, mpi_comm, 434191);

    auto n_cells = partitioned_grid_->boundaries[mpi_rank + 1]
                   - partitioned_grid_->boundaries[mpi_rank];

    auto gatherer = std::make_shared<AllVariablesGatherer>(
        std::move(cvars_gatherer), std::move(avars_gatherer), n_owned_cells);

    return std::make_shared<GatheredVisualization>(gatherer, single_node_vis);
  }

  std::shared_ptr<CFLCondition> choose_cfl_condition() override {
    auto cfl = super::choose_cfl_condition();

    auto op = ReductionOperation::min;
    auto all_reduce = std::make_shared<MPIAllreduce>(op, mpi_comm);

    return std::make_shared<DistributedCFLCondition>(cfl, all_reduce);
  }

  std::shared_ptr<BoundaryCondition> choose_boundary_condition() {




  }

protected:
  MPI_Comm mpi_comm;
  int mpi_rank;
  int mpi_comm_size;

  std::shared_ptr<Grid> full_grid_;
  std::shared_ptr<PartitionedGrid> partitioned_grid_;
};
}

#endif // ZISA_MPI_NUMERICAL_EXPERIMENT_HPP
