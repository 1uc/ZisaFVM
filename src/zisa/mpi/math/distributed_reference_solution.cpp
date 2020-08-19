#include <zisa/mpi/math/distributed_reference_solution.hpp>

#include <zisa/grid/neighbour_range.hpp>
#include <zisa/loops/reduction/all.hpp>
#include <zisa/math/max_quadrature_degree.hpp>
#include <zisa/mpi/io/hdf5_unstructured_writer.hpp>
#include <zisa/parallelization/distributed_grid.hpp>
#include <zisa/parallelization/domain_decomposition.hpp>

namespace zisa {

struct Sphere {
  XYZ x;
  double r;

  bool is_inside(const XYZ &p) { return zisa::norm(p - x) <= r; }
};

Sphere bounding_sphere(const Grid &grid) {
  auto x = XYZ::zeros();

  int_t count = 0;
  for (auto i : cell_indices(grid)) {
    if (grid.cell_flags[i].interior) {
      x += grid.cell_centers[i];
      count += 1;
    }
  }
  x = XYZ(x / double(count));

  double r = 0.0;
  for (auto i : cell_indices(grid)) {
    if (grid.cell_flags[i].interior) {
      for (auto k : neighbour_index_range(grid)) {
        r = max(r, zisa::norm(grid.vertex(i, k) - x));
      }
    }
  }

  return Sphere{x, r};
}

#define ZISA_TO_STRING(value, str)                                             \
  if (type == (value)) {                                                       \
    return (str);                                                              \
  }

std::string str(const DistributedReferenceSolution::OrgMessageType &type) {
  using E = DistributedReferenceSolution::OrgMessageType;

  ZISA_TO_STRING(E::ALL_DONE, "ALL_DONE");
  ZISA_TO_STRING(E::PE_IS_DONE, "PE_IS_DONE");
  ZISA_TO_STRING(E::XFER, "XFER");
  ZISA_TO_STRING(E::XFER_GOT_SOMETHING, "XFER_GOT_SOMETHING");

  LOG_ERR("Missing a case.");
}

#undef ZISA_TO_STRING

std::string str(const DistributedReferenceSolution::OrgMessage &msg) {
  return string_format("%s %d", str(msg.type).c_str(), msg.size);
}

DistributedReferenceSolution::DistributedReferenceSolution(
    serialize_t serialize,
    std::shared_ptr<Grid> large_grid,
    std::function<double(int_t, const XYZ &, int_t)> interpolation,
    int_t n_vars)
    : serialize(std::move(serialize)),
      large_grid(std::move(large_grid)),
      interpolation(std::move(interpolation)),
      n_vars(n_vars) {}

void DistributedReferenceSolution::compute_and_save(
    const std::string &output_name,
    const std::string &small_grid_folder,
    int new_small_comm_size,
    const std::function<bool(const Grid &, int_t)> &boundary_mask) {

  this->clear();

  std::vector<Sphere> bounding_spheres(integer_cast<size_t>(comm_size));
  bounding_spheres[integer_cast<size_t>(mpi_rank)]
      = bounding_sphere(*large_grid);
  zisa::mpi::allgather(array_view(bounding_spheres), comm);

  set_small_comm(new_small_comm_size);
  if (mpi_rank < small_comm_size) {
    auto small_grid_name = string_format("%s/%d/subgrid-%04d.msh.h5",
                                         small_grid_folder.c_str(),
                                         small_comm_size,
                                         mpi_rank);

    this->load_subgrid(small_grid_name, boundary_mask);

    auto n_cells = small_grid->n_cells;
    auto n_qp = small_grid->cells[0].qr.points.size();

    data = array<double, 3>(shape_t<3>{n_cells, n_vars, n_qp});

    completed_cells = array<bool, 2>(shape_t<2>{n_cells, n_qp});
    zisa::fill(completed_cells, false);

    for (auto [i, cell] : cells(*small_grid)) {
      if (small_grid->cell_flags[i].interior) {
        for (int_t k : index_range(n_qp)) {
          const auto &x = cell.qr.points[k];
          for (int p = 0; p < comm_size; ++p) {
            if (bounding_spheres[size_t(p)].is_inside(x)) {
              ids[p].push_back({i, k});
              coords[p].push_back(x);
            }
          }
        }
      } else {
        for (int_t k : index_range(n_qp)) {
          for (int_t l : index_range(n_vars)) {
            data(i, l, k) = 0.0;
          }
          completed_cells(i, k) = true;
        }
      }
    }

    for (const auto &[p, id] : ids) {
      const auto &coord = coords[p];
      send_org_message(
          OrgMessageType::XFER, integer_cast<int>(coord.size()), p);

      send_xfer_request(coord, p);
    }
  }

  while (true) {
    auto [msg, src] = receive_org_message();

    if (msg.type == OrgMessageType::ALL_DONE) {
      break;
    } else if (msg.type == OrgMessageType::PE_IS_DONE) {
      ++n_pes_done;
      if (all_done()) {
        for (int dst = 0; dst < comm_size; ++dst) {
          send_org_message(OrgMessageType::ALL_DONE, 0, dst);
        }
      }

    } else if (msg.type == OrgMessageType::XFER) {
      process_xfer_request(msg, src);
    } else if (msg.type == OrgMessageType::XFER_GOT_SOMETHING) {
      process_xfer_response(msg, src);

      if (is_done()) {
        send_org_message(
            OrgMessageType::PE_IS_DONE, /* msg_size */ 0, /* dst */ 0);
      }
    }
  }

  zisa::mpi::wait_all(requests);
  zisa::mpi::barrier(comm);

  if (mpi_rank < small_comm_size) {
    auto all_vars = AllVariables{};
    all_vars.cvars = average_data();

    std::string stem = zisa::basename(small_grid_folder);
    std::string filename = string_format(
        "down_sampled/%s/%s", stem.c_str(), output_name.c_str());

    serialize(filename, small_comm, *small_dgrid, all_vars);
  }
}

void DistributedReferenceSolution::set_small_comm(int new_small_comm_size) {
  MPI_Comm_split(comm, mpi_rank < new_small_comm_size, mpi_rank, &small_comm);
  this->small_comm_size = new_small_comm_size;
}

void DistributedReferenceSolution::load_subgrid(
    const std::string &grid_name,
    const std::function<bool(const Grid &, int_t)> &mask) {
  auto grid = load_grid(grid_name);
  auto dgrid = load_distributed_grid(grid_name);

  auto is_inside = [this, &dgrid](int_t i) {
    return dgrid.partition[i] == integer_cast<int_t>(mpi_rank);
  };

  auto [vertex_indices, vertices, global_indices]
      = extract_subgrid_v2(*grid, is_inside);

  small_grid = std::make_shared<Grid>(grid->element_type(),
                                      std::move(vertices),
                                      std::move(vertex_indices),
                                      MAX_QUADRATURE_DEGREE);

  mask_ghost_cells(*small_grid, mask);

  small_dgrid = std::make_shared<DistributedGrid>(
      extract_distributed_subgrid(dgrid, global_indices));
}

double DistributedReferenceSolution::interpolate_reference_solution(
    int_t i, const XYZ &x, int_t k) {
  return interpolation(i, x, k);
}

GridVariables DistributedReferenceSolution::average_data() {
  auto cvars = GridVariables(shape_t<2>{small_grid->n_cells, n_vars});
  for (const auto &[i, cell] : cells(*small_grid)) {
    const auto qr = cell.qr;
    auto n_qr = qr.points.size();
    cvars(i, 0) = 0.0;
    for (int_t l : index_range(n_vars)) {
      for (int_t k : index_range(n_qr)) {
        cvars(i, l) += qr.weights[k] / qr.volume * data(i, l, k);
      }
    }
  }

  return cvars;
}

bool DistributedReferenceSolution::all_done() {
  return n_pes_done == small_comm_size;
}

bool DistributedReferenceSolution::is_done() {
  return zisa::reduce::all(flat_range(completed_cells),
                           [this](int_t i) { return completed_cells[i]; });
}

std::pair<DistributedReferenceSolution::OrgMessage, int>
DistributedReferenceSolution::receive_org_message() {

  auto [msg, status]
      = zisa::mpi::recv_pod<OrgMessage>(MPI_ANY_SOURCE, org_tag, comm);

  return {msg, status.source};
}

void DistributedReferenceSolution::process_xfer_request(
    const DistributedReferenceSolution::OrgMessage &msg, int dst) {
  auto coordinates = zisa::array<XYZ, 1>(msg.size);
  zisa::mpi::recv(array_view(coordinates), dst, xfer_tag_coords, comm);

  auto indices = zisa::array<int_t, 1>(msg.size);
  auto values = zisa::array<double, 2>(shape_t<2>{msg.size, n_vars});

  int_t count = 0;
  for (int_t k : index_range(values.shape(0))) {
    const XYZ &x = coordinates[k];
    auto o = locate(*large_grid, x);
    if (o == std::nullopt) {
      continue;
    }

    int_t i = o.value();
    if (large_grid->cell_flags[i].interior) {
      indices[count] = k;
      for (int_t var : index_range(n_vars)) {
        values(count, var) = interpolate_reference_solution(i, x, var);
      }
      ++count;
    }
  }

  // For `count == 0` nothing needs to be done.
  if (count > 0) {
    indices_vault.push_back(std::move(indices));
    values_vault.push_back(std::move(values));

    send_org_message(
        OrgMessageType::XFER_GOT_SOMETHING, integer_cast<int>(count), dst);
    auto indices_request = zisa::mpi::isend(
        const_slice(array_const_view(indices_vault.back()), 0, count),
        dst,
        xfer_tag_indices,
        comm);
    requests.push_back(std::move(indices_request));

    auto values_request = zisa::mpi::isend(
        const_slice(array_const_view(values_vault.back()), 0, count),
        dst,
        xfer_tag_values,
        comm);
    requests.push_back(std::move(values_request));
  }
}

void DistributedReferenceSolution::process_xfer_response(
    const DistributedReferenceSolution::OrgMessage &msg, int src) {
  auto indices = array<int_t, 1>(msg.size);
  auto values = array<double, 2>(shape_t<2>{msg.size, n_vars});

  zisa::mpi::recv(array_view(indices), src, xfer_tag_indices, comm);
  zisa::mpi::recv(array_view(values), src, xfer_tag_values, comm);

  const auto &ids = this->ids.at(src);
  for (int_t k = 0; k < indices.size(); ++k) {
    auto [i, l] = ids.at(indices[k]);

    for (int_t var : index_range(n_vars)) {
      data(i, var, l) = values(k, var); // make multi dim.
    }

    completed_cells(i, l) = true;
  }
}

void DistributedReferenceSolution::send_org_message(
    const DistributedReferenceSolution::OrgMessageType &msg_type,
    int msg_size,
    int dst) {

  org_messages_vault.emplace_back(msg_type, msg_size);

  requests.push_back(
      zisa::mpi::isend_pod(org_messages_vault.back(), dst, org_tag, comm));
}

void DistributedReferenceSolution::send_xfer_request(
    const std::vector<XYZ> &coord, int dst) {
  requests.push_back(
      zisa::mpi::isend(array_const_view(coord), dst, xfer_tag_coords, comm));
}

void DistributedReferenceSolution::clear() {
  requests.clear();
  org_messages_vault.clear();
  indices_vault.clear();
  values_vault.clear();
  ids.clear();
  coords.clear();
  n_pes_done = 0;

  if (small_comm_size != -1) {
    MPI_Comm_free(&small_comm);
    small_comm_size = -1;
  }
}
}