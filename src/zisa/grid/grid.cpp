#include <algorithm>
#include <map>
#include <vector>

#include <zisa/config.hpp>
#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/grid/grid_impl.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/utils/logging.hpp>

namespace zisa {
int_t count_edges(const neighbours_t &neighbours, const is_valid_t &is_valid) {
  int_t n_edges = 0;
  int_t n_cells = neighbours.shape(0);
  int_t max_neighbours = neighbours.shape(1);

  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {

      int_t j = neighbours(i, k);
      if (!is_valid(i, k) || i > j) {
        continue;
      }

      ++n_edges;
    }
  }

  return n_edges;
}

normals_t compute_normals(const vertices_t &vertices,
                          const vertex_indices_t &vertex_indices,
                          const neighbours_t &neighbours,
                          const is_valid_t &is_valid,
                          const edge_indices_t &edge_indices) {

  auto n_edges = count_edges(neighbours, is_valid);
  auto n_cells = vertices.shape(0);
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
    const auto &v0 = vertices[vertex_indices(i, 0)];
    const auto &v1 = vertices[vertex_indices(i, 1)];
    const auto &v2 = vertices[vertex_indices(i, 2)];

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

  auto cell_centers = empty_like(vertices);

  int_t n_cells = vertices.shape(0);
  for (int_t i = 0; i < n_cells; ++i) {

    const auto &v1 = vertices[vertex_indices(i, 0)];
    const auto &v2 = vertices[vertex_indices(i, 1)];
    const auto &v3 = vertices[vertex_indices(i, 2)];

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

  int_t n_edges = 0;
  auto n_cells = neighbours.shape(0);
  auto max_neighbours = neighbours.shape(1);
  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {

      if (!is_valid(i, k)) {
        edge_indices(i, k) = magic_index_value;
      } else {
        auto j = neighbours(i, k);
        if (i < j) {
          edge_indices(i, k) = n_edges++;
        } else {
          auto kj = find_self(neighbours, i, j);
          edge_indices(i, k) = edge_indices(j, kj);
        }
      }
    }
  }

  return edge_indices;
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

edges_t compute_edges(const vertex_indices_t &vertex_indices) {
  auto n_vertices = count_vertices(vertex_indices);
  auto n_cells = vertex_indices.shape(0);

  edges_t ret(n_vertices);

  assert(vertex_indices.shape(0) == n_cells);
  for (int_t i = 0; i < n_cells; ++i) {

    int_t v0 = vertex_indices(i, 0);
    int_t v1 = vertex_indices(i, 1);
    int_t v2 = vertex_indices(i, 2);

    ret[v0][v1] = i;
    ret[v1][v2] = i;
    ret[v2][v0] = i;
  }

  return ret;
}

namespace /* anon */ {
int_t first_empty(const neighbours_t &neighbours, int_t i) {
  int_t k = 0;
  while (neighbours(i, k) != magic_index_value) {
    ++k;
  }
  return k;
};
} // namespace

neighbours_t compute_neighbours(const vertex_indices_t &vertex_indices) {
  // Only works when all polygons in the grid have the same number of edges.
  assert(vertex_indices.shape(1) == 3);

  int_t n_cells = vertex_indices.shape(0);
  int_t max_neighbours = vertex_indices.shape(1);

  auto neighbours = empty_like(vertex_indices);
  std::fill(neighbours.begin(), neighbours.end(), magic_index_value);

  auto edges = compute_edges(vertex_indices);

  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {
      const auto &v0 = vertex_indices(i, k);
      const auto &v1 = vertex_indices(i, (k + 1) % max_neighbours);

      if (edges[v1].find(v0) != edges[v1].end()) {
        neighbours(i, k) = edges[v1][v0];
      }
    }
  }

  return neighbours;
}

Grid::Grid(array<XY, 1> vertices_, array<int_t, 2> vertex_indices_)
    : vertices(std::move(vertices_)), vertex_indices(std::move(vertex_indices_)) {

  n_cells = vertex_indices.shape(0);
  max_neighbours = vertex_indices.shape(1);

  neighbours = compute_neighbours(this->vertex_indices);
  is_valid = compute_valid_neighbours(neighbours);

  cell_centers = compute_cell_centers(this->vertices, this->vertex_indices);
  edge_indices = compute_edge_indices(neighbours, is_valid);
  volumes = compute_volumes(this->vertices, this->vertex_indices);

  normals = compute_normals(
      this->vertices, this->vertex_indices, neighbours, is_valid, edge_indices);

  tangentials = compute_tangentials(normals);
}

Triangle Grid::triangles(int_t i) const {
  const auto &v0 = vertices[vertex_indices(i, 0)];
  const auto &v1 = vertices[vertex_indices(i, 1)];
  const auto &v2 = vertices[vertex_indices(i, 2)];

  return Triangle(v0, v1, v2);
}

std::shared_ptr<Grid> load_gmsh(const std::string &filename) {

  auto gmsh = GMSHReader(filename);

  auto max_neighbours = int_t(3);
  auto n_vertices = gmsh.vertices.size();
  auto n_cells = gmsh.vertex_indices.size();

  auto vertices = array<XY, 1>(shape_t<1>{n_vertices});
  auto vertex_indices = array<int_t, 2>(shape_t<2>{n_cells, max_neighbours});

  for (int_t i = 0; i < n_vertices; ++i) {
    const auto &v = gmsh.vertices[i];
    vertices(i) = {v[0], v[1]};
  }

  for (int_t i = 0; i < n_cells; ++i) {
    const auto &v = gmsh.vertex_indices[i];
    for (int_t k = 0; k < max_neighbours; ++k) {
      vertex_indices(i, k) = v[k];
    }
  }

  return std::make_shared<Grid>(std::move(vertices), std::move(vertex_indices));
}

double largest_circum_radius(const Grid &grid) {

  double r = 0.0;
  for (auto &&tri : triangles(grid)) {
    r = zisa::max(r, circum_radius(tri));
  }

  return r;
}

} // namespace zisa
