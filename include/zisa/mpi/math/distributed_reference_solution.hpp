// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_DISTRIBUTED_REFERENCE_SOLUTION_HPP_OQYEB
#define ZISA_DISTRIBUTED_REFERENCE_SOLUTION_HPP_OQYEB

#include <zisa/config.hpp>

#include <functional>
#include <list>
#include <zisa/grid/grid.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/mpi/io/hdf5_unstructured_writer.hpp>
#include <zisa/mpi/mpi.hpp>
#include <zisa/mpi/mpi_tag_constants.hpp>
#include <zisa/parallelization/distributed_grid.hpp>

namespace zisa {

class DistributedReferenceSolution {
public:
  enum class OrgMessageType { XFER, XFER_GOT_SOMETHING, PE_IS_DONE, ALL_DONE };

  struct OrgMessage {
    OrgMessageType type;
    int size;

    OrgMessage() = default;
    OrgMessage(const OrgMessage &) = default;
    OrgMessage(const OrgMessageType &msg, int size) : type(msg), size(size) {}
  };

private:
  using serialize_t = std::function<void(
      std::string, MPI_Comm, const DistributedGrid &, const AllVariables)>;

public:
  DistributedReferenceSolution(
      serialize_t serialize,
      std::shared_ptr<Grid> large_grid,
      std::function<double(int_t, const XYZ &, int_t)> interpolation,
      int_t n_vars);

  void compute_and_save(
      const std::string &output_name,
      const std::string &small_grid_folder,
      int small_comm_size,
      const std::function<bool(const Grid &grid, int_t)> &boundary_mask);

protected:
  void clear();
  void set_small_comm(int small_comm_size);

  void load_subgrid(const std::string &grid_name,
                    const std::function<bool(const Grid &grid, int_t)> &mask);

  double interpolate_reference_solution(int_t i, const XYZ &x, int_t k);

  GridVariables average_data();

  bool all_done();
  bool is_done();

  void send_xfer_request(const std::vector<XYZ> &coord, int dst);
  void process_xfer_request(const OrgMessage &msg, int dst);
  void process_xfer_response(const OrgMessage &msg, int src);

  void send_org_message(const OrgMessageType &msg, int size, int dst);
  std::pair<OrgMessage, int> receive_org_message();

private:
  serialize_t serialize;
  std::shared_ptr<Grid> large_grid;
  std::function<double(int_t, const XYZ &, int_t)> interpolation;
  int_t n_vars;

  MPI_Comm comm = MPI_COMM_WORLD;
  int comm_size = zisa::mpi::size();
  int mpi_rank = zisa::mpi::rank();

  MPI_Comm small_comm;
  int small_comm_size = -1;

  int xfer_tag_coords = ZISA_MPI_TAG_REFERENCE_XFER;
  int xfer_tag_indices = ZISA_MPI_TAG_REFERENCE_XFER + 1;
  int xfer_tag_values = ZISA_MPI_TAG_REFERENCE_XFER + 2;
  int org_tag = ZISA_MPI_TAG_REFERENCE_ORG;

  std::vector<zisa::mpi::Request> requests;

  std::list<OrgMessage> org_messages_vault;
  std::list<array<int_t, 1>> indices_vault;
  std::list<array<double, 2>> values_vault;
  std::unordered_map<int, std::vector<std::pair<int_t, int_t>>> ids;
  std::unordered_map<int, std::vector<XYZ>> coords;

  std::shared_ptr<Grid> small_grid;
  std::shared_ptr<DistributedGrid> small_dgrid;

  array<double, 3> data;
  array<bool, 2> completed_cells;

  int n_pes_done = 0;
};

}

#endif // ZISA_DISTRIBUTED_REFERENCE_SOLUTION_HPP
