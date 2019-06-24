#include <zisa/grid/grid.hpp>

#include <algorithm>
#include <map>
#include <optional>
#include <vector>

#include <zisa/config.hpp>
#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/grid/grid_impl.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/io/hdf5_writer.hpp>
#include <zisa/loops/reduction/max.hpp>
#include <zisa/loops/reduction/min.hpp>
#include <zisa/loops/reduction/sum.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/math/symmetric_choices.hpp>
#include <zisa/utils/logging.hpp>

namespace zisa {
int_t count_interior_edges(const neighbours_t &neighbours,
                           const is_valid_t &is_valid) {
  int_t n_interior_edges = 0;
  int_t n_cells = neighbours.shape(0);
  int_t max_neighbours = neighbours.shape(1);

  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {

      int_t j = neighbours(i, k);
      if (!is_valid(i, k) || i > j) {
        continue;
      }

      ++n_interior_edges;
    }
  }

  return n_interior_edges;
}

int_t count_exterior_edges(const is_valid_t &is_valid) {
  auto n_cells = is_valid.shape(0);
  auto max_neighbours = is_valid.shape(1);

  int_t n_exterior_edges = 0;
  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {
      if (!is_valid(i, k)) {
        ++n_exterior_edges;
      }
    }
  }

  return n_exterior_edges;
}

normals_t compute_normals(const vertices_t &vertices,
                          const vertex_indices_t &vertex_indices,
                          const neighbours_t &neighbours,
                          const is_valid_t &is_valid,
                          const edge_indices_t &edge_indices) {

  auto n_interior_edges = count_interior_edges(neighbours, is_valid);
  auto n_exterior_edges = count_exterior_edges(is_valid);
  auto n_edges = n_interior_edges + n_exterior_edges;

  auto n_cells = vertex_indices.shape(0);
  auto max_neighbours = vertex_indices.shape(1);

  auto normals = normals_t(n_edges);
  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {

      int_t j = neighbours(i, k);
      if (!is_valid(i, k) || i > j) {
        continue;
      }

      const auto &vi = vertices[vertex_indices(i, k)];
      const auto &vj = vertices[vertex_indices(i, (k + 1) % max_neighbours)];

      int_t ei = edge_indices(i, k);
      normals(ei) = rotate_right(normalize(vj - vi));
    }
  }

  return normals;
}

tangentials_t compute_tangentials(const normals_t &normals) {
  auto tangentials = empty_like(normals);

  int_t n_edges = tangentials.shape(0);
  for (int_t ei = 0; ei < n_edges; ++ei) {
    tangentials(ei) = rotate_left(normals(ei));
  }

  return tangentials;
}

volumes_t compute_volumes(const vertices_t &vertices,
                          const vertex_indices_t &vertex_indices) {
  assert(vertex_indices.shape(1) == 3);

  int_t n_cells = vertex_indices.shape(0);
  auto volumes = volumes_t(shape_t<1>(n_cells));

  for (int_t i = 0; i < n_cells; ++i) {
    const auto &v0 = vertices[vertex_indices(i, int_t(0))];
    const auto &v1 = vertices[vertex_indices(i, int_t(1))];
    const auto &v2 = vertices[vertex_indices(i, int_t(2))];

    auto a = norm(v1 - v0);
    auto b = norm(v2 - v1);
    auto c = norm(v0 - v2);

    volumes(i) = herons_formula(a, b, c);
  }

  return volumes;
}

cell_centers_t compute_cell_centers(const vertices_t &vertices,
                                    const vertex_indices_t &vertex_indices) {
  assert(vertex_indices.shape(1) == 3);

  int_t n_cells = vertex_indices.shape(0);
  auto cell_centers = cell_centers_t(shape_t<1>{n_cells});
  for (int_t i = 0; i < n_cells; ++i) {

    const auto &v1 = vertices[vertex_indices(i, int_t(0))];
    const auto &v2 = vertices[vertex_indices(i, int_t(1))];
    const auto &v3 = vertices[vertex_indices(i, int_t(2))];

    cell_centers(i) = (v1 + v2 + v3) / 3.0;
  }

  return cell_centers;
}

int_t find_self(const neighbours_t &neighbours, int_t me, int_t other) {

  auto max_neighbours = neighbours.shape(1);
  for (int_t k = 0; k < max_neighbours; ++k) {
    if (neighbours(other, k) == me) {
      return k;
    }
  }

  LOG_ERR("Failed to find myself.");
}

edge_indices_t compute_edge_indices(const neighbours_t &neighbours,
                                    const is_valid_t &is_valid) {

  auto edge_indices = empty_like(neighbours);

  int_t n_interior_edges = 0;
  int_t n_exterior_edges = count_interior_edges(neighbours, is_valid);

  auto n_cells = neighbours.shape(0);
  auto max_neighbours = neighbours.shape(1);

  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {

      if (!is_valid(i, k)) {
        edge_indices(i, k) = n_exterior_edges++;
      } else {
        auto j = neighbours(i, k);
        if (i < j) {
          edge_indices(i, k) = n_interior_edges++;
        } else {
          auto kj = find_self(neighbours, i, j);
          edge_indices(i, k) = edge_indices(j, kj);
        }
      }
    }
  }

  return edge_indices;
}

left_right_t compute_left_right(const edge_indices_t &edge_indices,
                                const neighbours_t &neighbours,
                                const is_valid_t &is_valid) {
  int_t n_cells = edge_indices.shape(0);
  int_t max_neighbours = edge_indices.shape(1);
  int_t n_interior_edges = count_interior_edges(neighbours, is_valid);
  int_t n_exterior_edges = count_exterior_edges(is_valid);
  int_t n__edges = n_interior_edges + n_exterior_edges;

  auto left_right = left_right_t(shape_t<1>{n__edges});

  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {
      auto e = edge_indices(i, k);

      auto j = neighbours(i, k);
      if (is_valid(i, k)) {
        if (i < j) {
          left_right(e) = std::pair<int_t, int_t>(i, j);
        }
      } else {
        left_right(e) = std::pair<int_t, int_t>(i, magic_index_value);
      }
    }
  }

  return left_right;
}

is_valid_t compute_valid_neighbours(const neighbours_t &neighbours) {
  auto is_valid = array<bool, 2>(neighbours.shape());

  auto n_cells = neighbours.shape(0);
  auto max_neighbours = neighbours.shape(1);
  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {
      is_valid(i, k) = (neighbours(i, k) != magic_index_value);
    }
  }

  return is_valid;
}

int_t count_vertices(const vertex_indices_t &vertex_indices) {
  return *std::max_element(vertex_indices.begin(), vertex_indices.end()) + 1;
}

vertex_neighbours_t
compute_vertex_neighbours(const vertex_indices_t &vertex_indices) {
  auto n_vertices = count_vertices(vertex_indices);
  auto n_cells = vertex_indices.shape(0);
  auto max_neighbours = vertex_indices.shape(1);

  vertex_neighbours_t ret(n_vertices);

  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {

      auto v = vertex_indices(i, k);
      ret[v].push_back(i);
    }
  }

  return ret;
}

std::optional<std::pair<int_t, int_t>>
common_face(const vertex_indices_t &vertex_indices, int_t i, int_t j) {
  // TODO this only works for triangles.

  auto max_neighbours = vertex_indices.shape(1);
  auto si = std::vector<bool>(max_neighbours);
  auto sj = std::vector<bool>(max_neighbours);

  std::fill(si.begin(), si.end(), false);
  std::fill(sj.begin(), sj.end(), false);

  for (int_t k = 0; k < max_neighbours; ++k) {
    auto vi = vertex_indices(i, k);

    for (int_t l = 0; l < max_neighbours; ++l) {
      if (vi == vertex_indices(j, l)) {
        si[k] = true;
        sj[l] = true;
      }
    }
  }

  // This is the part that is specific for triangles.
  auto relative_neighbour_index = [](const auto &s) {
    for (int_t k = 0; k < 3; ++k) {
      if (s[k] && s[(k + 1) % 3]) {
        return k;
      }
    }

    return magic_index_value;
  };

  auto ki = relative_neighbour_index(si);
  auto kj = relative_neighbour_index(sj);

  if (ki == magic_index_value) {
    return std::nullopt;
  }

  return std::optional<std::pair<int_t, int_t>>({ki, kj});
}

neighbours_t compute_neighbours(const vertex_indices_t &vertex_indices) {
  int_t n_vertices = count_vertices(vertex_indices);

  auto neighbours = empty_like(vertex_indices);
  std::fill(neighbours.begin(), neighbours.end(), magic_index_value);

  auto vertex_neighbours = compute_vertex_neighbours(vertex_indices);

  for (int_t vi = 0; vi < n_vertices; ++vi) {
    for (auto [i, j] : strict_symmetric_choices(vertex_neighbours[vi])) {
      auto face = common_face(vertex_indices, i, j);
      if (face) {
        neighbours(i, face->first) = j;
        neighbours(j, face->second) = i;
      }
    }
  }

  return neighbours;
}

array<array<double, 1>, 1> compute_normalized_moments(const Grid &grid) {
  auto n_cells = grid.n_cells;
  auto normalized_moments = array<array<double, 1>, 1>(shape_t<1>{n_cells});
  for (const auto &[i, tri] : triangles(grid)) {
    normalized_moments(i) = zisa::normalized_moments(tri, 4, 4);
  }

  return normalized_moments;
}

Grid::Grid(array<XYZ, 1> vertices_, array<int_t, 2> vertex_indices_)
    : vertex_indices(std::move(vertex_indices_)),
      vertices(std::move(vertices_)) {

  n_cells = vertex_indices.shape(0);
  n_vertices = vertices.shape(0);
  max_neighbours = vertex_indices.shape(1);

  neighbours = compute_neighbours(this->vertex_indices);
  is_valid = compute_valid_neighbours(neighbours);

  n_interior_edges = count_interior_edges(neighbours, is_valid);
  n_exterior_edges = count_exterior_edges(is_valid);
  n_edges = n_interior_edges + n_exterior_edges;

  cell_centers = compute_cell_centers(this->vertices, this->vertex_indices);
  edge_indices = compute_edge_indices(neighbours, is_valid);
  left_right = compute_left_right(edge_indices, neighbours, is_valid);

  volumes = compute_volumes(this->vertices, this->vertex_indices);

  normals = compute_normals(
      this->vertices, this->vertex_indices, neighbours, is_valid, edge_indices);

  tangentials = compute_tangentials(normals);

  normalized_moments = compute_normalized_moments(*this);
}

const XYZ &Grid::vertex(int_t i, int_t k) const {
  return vertices(vertex_indices(i, k));
}

Triangle Grid::triangle(int_t i) const {
  const auto &v0 = vertex(i, 0);
  const auto &v1 = vertex(i, 1);
  const auto &v2 = vertex(i, 2);

  return Triangle(v0, v1, v2);
}

Edge Grid::edge(int_t e) const {
  int_t i = left_right(e).first;

  for (int_t k = 0; k < max_neighbours - 1; ++k) {
    if (edge_indices(i, k) == e) {
      return Edge(vertex(i, k), vertex(i, k + 1));
    }
  }

  return Edge(vertex(i, max_neighbours - 1), vertex(i, int_t(0)));
}

Edge Grid::edge(int_t i, int_t k) const { return edge(edge_indices(i, k)); }

double Grid::characteristic_length(int_t i) const {
  return zisa::characteristic_length(triangle(i));
}

std::string Grid::str() const {
  double dx_min = smallest_inradius(*this);
  double dx_max = largest_circum_radius(*this);
  return string_format("n_cells : %d\n"
                       "n_vertices : %d\n"
                       "n_edges : %d\n"
                       "dx_min : %e\n"
                       "dx_max : %e\n",
                       n_cells,
                       n_vertices,
                       n_edges,
                       dx_min,
                       dx_max);
}

std::optional<int_t>
locate(const Grid &grid, const XYZ &x, int_t i_guess, int_t max_iter) {
  auto max_neighbours = grid.max_neighbours;

  int_t i = i_guess;

  for (int_t count = 0; count < max_iter; ++count) {
    if (is_inside(grid.triangle(i), x)) {
      return i;
    } else {

      auto connecting_edge = Edge(grid.cell_centers(i), x);
      for (int_t k = 0; k < max_neighbours; ++k) {
        if (grid.is_valid(i, k)
            && is_intersecting(grid.edge(i, k), connecting_edge)) {

          i = grid.neighbours(i, k);
          break;
        }
      }
    }
  }

  return std::nullopt;
}

std::optional<int_t> locate(const Grid &grid, const XYZ &x) {
  return locate(grid, x, 0, grid.n_cells);
}

double volume(const Grid &grid) {
  return zisa::reduce::sum(
      triangles(grid), [](int_t, const Triangle &tri) { return volume(tri); });
}

std::shared_ptr<Grid> load_grid(const std::string &filename) {

  auto len = filename.size();

  if (filename.substr(len - 4) == ".msh") {
    return load_gmsh(filename);
  } else if (filename.substr(len - 3) == ".h5") {
    auto reader = HDF5SerialReader(filename);
    return std::make_shared<Grid>(Grid::load(reader));
  }

  LOG_ERR(string_format("Unknown filetype. [%s]", filename.c_str()));
}

std::shared_ptr<Grid> load_gmsh(const std::string &filename) {
  auto gmsh = GMSHData(filename);

  auto max_neighbours = int_t(3);
  auto n_vertices = gmsh.vertices.size();
  auto n_cells = gmsh.vertex_indices.size();

  auto vertices = array<XYZ, 1>(shape_t<1>{n_vertices});
  auto vertex_indices = array<int_t, 2>(shape_t<2>{n_cells, max_neighbours});

  for (int_t i = 0; i < n_vertices; ++i) {
    const auto &v = gmsh.vertices[i];
    vertices(i) = {v[0], v[1], v[2]};
  }

  for (int_t i = 0; i < n_cells; ++i) {
    const auto &v = gmsh.vertex_indices[i];
    for (int_t k = 0; k < max_neighbours; ++k) {
      vertex_indices(i, k) = v[k];
    }
  }

  return std::make_shared<Grid>(std::move(vertices), std::move(vertex_indices));
}

void save(HDF5Writer &writer, const Grid &grid) {
  writer.write_scalar(grid.n_cells, "n_cells");
  writer.write_scalar(grid.n_vertices, "n_vertices");
  writer.write_scalar(grid.n_edges, "n_edges");
  writer.write_scalar(grid.n_interior_edges, "n_interior_edges");
  writer.write_scalar(grid.n_exterior_edges, "n_exterior_edges");
  writer.write_scalar(grid.max_neighbours, "max_neighbours");

  save(writer, grid.vertex_indices, "vertex_indices");
  save(writer, grid.edge_indices, "edge_indices");

  save(writer, grid.neighbours, "neighbours");
  save(writer, grid.is_valid, "is_valid");

  save(writer, grid.vertices, "vertices");
  save(writer, grid.cell_centers, "cell_centers");

  save(writer, grid.volumes, "volumes");
  save(writer, grid.normals, "normals");
  save(writer, grid.tangentials, "tangentials");

  writer.write_scalar(largest_circum_radius(grid), "dx_max");
  writer.write_scalar(smallest_inradius(grid), "dx_min");
}

double largest_circum_radius(const Grid &grid) {
  return zisa::reduce::max(triangles(grid), [](int_t, const Triangle &tri) {
    return circum_radius(tri);
  });
}

double smallest_inradius(const Grid &grid) {
  return zisa::reduce::min(triangles(grid), [](int_t, const Triangle &tri) {
    return inradius(tri);
  });
}

Grid Grid::load(HDF5Reader &reader) {
  auto grid = Grid{};

  grid.n_cells = reader.read_scalar<int_t>("n_cells");
  grid.n_vertices = reader.read_scalar<int_t>("n_vertices");
  grid.n_edges = reader.read_scalar<int_t>("n_edges");
  grid.n_interior_edges = reader.read_scalar<int_t>("n_interior_edges");
  grid.n_exterior_edges = reader.read_scalar<int_t>("n_exterior_edges");
  grid.max_neighbours = reader.read_scalar<int_t>("max_neighbours");

  grid.vertex_indices = array<int_t, 2>::load(reader, "vertex_indices");
  grid.edge_indices = array<int_t, 2>::load(reader, "edge_indices");

  grid.neighbours = array<int_t, 2>::load(reader, "neighbours");
  grid.is_valid = array<bool, 2>::load(reader, "is_valid");

  grid.vertices = array<XYZ, 1>::load(reader, "vertices");
  grid.cell_centers = array<XYZ, 1>::load(reader, "cell_centers");

  grid.volumes = array<double, 1>::load(reader, "volumes");
  grid.normals = array<XYZ, 1>::load(reader, "normals");
  grid.tangentials = array<XYZ, 1>::load(reader, "tangentials");

  grid.left_right
      = compute_left_right(grid.edge_indices, grid.neighbours, grid.is_valid);
  grid.normalized_moments = compute_normalized_moments(grid);

  return grid;
}

array<double, 1>
normalized_moments(const Triangle &tri, int degree, int_t quad_deg) {
  auto m = array<double, 1>(shape_t<1>{poly_dof(degree)});

  auto length = characteristic_length(tri);
  double length_d = 1.0;

  for (int d = 0; d <= degree; ++d) {
    for (int k = 0; k <= d; ++k) {
      int l = d - k;

      m(poly_index(k, l)) = avg_moment(tri, k, l, quad_deg) / length_d;
    }

    length_d *= length;
  }

  return m;
}

} // namespace zisa
