#include <zisa/grid/grid.hpp>

#include <algorithm>
#include <list>
#include <map>
#include <optional>
#include <vector>

#include <zisa/config.hpp>
#include <zisa/grid/cell_range.hpp>
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
#include <zisa/math/tetrahedral_rule.hpp>
#include <zisa/math/tetrahedron.hpp>
#include <zisa/math/triangular_rule.hpp>
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

normals_t compute_normals(GMSHElementType element_type,
                          const vertices_t &vertices,
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

    auto vertex_index = [&vertex_indices, &element_type, i](int_t face,
                                                            int_t rel) {
      auto k = GMSHElementInfo::relative_vertex_index(element_type, face, rel);
      return vertex_indices(i, k);
    };

    for (int_t k = 0; k < max_neighbours; ++k) {

      int_t j = neighbours(i, k);
      if (!is_valid(i, k) || i > j) {
        continue;
      }

      int_t ei = edge_indices(i, k);
      if (element_type == GMSHElementType::triangle) {
        const auto &v0 = vertices[vertex_index(k, 0)];
        const auto &v1 = vertices[vertex_index(k, 1)];

        normals(ei) = rotate_right(normalize(v1 - v0));
      } else if (element_type == GMSHElementType::tetrahedron) {
        const auto &v0 = vertices[vertex_index(k, 0)];
        const auto &v1 = vertices[vertex_index(k, 1)];
        const auto &v2 = vertices[vertex_index(k, 2)];

        normals(ei) = XYZ(zisa::normalize(zisa::cross(v1 - v0, v2 - v0)));
      }
    }
  }

  return normals;
}

tangentials_t compute_tangentials(GMSHElementType element_type,
                                  const normals_t &normals,
                                  const vertices_t &vertices,
                                  const vertex_indices_t &vertex_indices,
                                  const neighbours_t &neighbours,
                                  const is_valid_t &is_valid,
                                  const edge_indices_t &edge_indices) {
  auto n_edges = normals.shape(0);
  auto n_cells = vertex_indices.shape(0);
  auto max_neighbours = vertex_indices.shape(1);

  auto tangentials = tangentials_t(shape_t<2>{n_edges, int_t(2)});
  for (int_t i = 0; i < n_cells; ++i) {

    auto vertex_index = [&vertex_indices, &element_type, i](int_t face,
                                                            int_t rel) {
      auto k = GMSHElementInfo::relative_vertex_index(element_type, face, rel);
      return vertex_indices(i, k);
    };

    for (int_t k = 0; k < max_neighbours; ++k) {

      int_t j = neighbours(i, k);
      if (!is_valid(i, k) || i > j) {
        continue;
      }

      int_t ei = edge_indices(i, k);
      const auto &v0 = vertices[vertex_index(k, int_t(0))];
      const auto &v1 = vertices[vertex_index(k, int_t(1))];

      tangentials(ei, int_t(0)) = zisa::normalize(v1 - v0);
      tangentials(ei, int_t(1))
          = zisa::cross(normals(ei), tangentials(ei, int_t(0)));
    }
  }

  return tangentials;
}

double volume(GMSHElementType element_type,
              const vertices_t &vertices,
              const vertex_indices_t &vertex_indices,
              int_t i) {

  const auto &v0 = vertices[vertex_indices(i, int_t(0))];
  const auto &v1 = vertices[vertex_indices(i, int_t(1))];
  const auto &v2 = vertices[vertex_indices(i, int_t(2))];

  if (element_type == GMSHElementType::triangle) {
    return volume(Triangle{v0, v1, v2});
  } else if (element_type == GMSHElementType::tetrahedron) {
    const auto &v3 = vertices[vertex_indices(i, int_t(3))];
    return volume(Tetrahedron{v0, v1, v2, v3});
  }

  LOG_ERR("Unknown element type.")
}

volumes_t compute_volumes(GMSHElementType element_type,
                          const vertices_t &vertices,
                          const vertex_indices_t &vertex_indices) {
  int_t n_cells = vertex_indices.shape(0);
  auto volumes = volumes_t(shape_t<1>(n_cells));

  for (int_t i = 0; i < n_cells; ++i) {
    volumes(i) = volume(element_type, vertices, vertex_indices, i);
  }

  return volumes;
}

XYZ cell_center(GMSHElementType element_type,
                const vertices_t &vertices,
                const vertex_indices_t &vertex_indices,
                int_t i) {
  const auto &v0 = vertices[vertex_indices(i, int_t(0))];
  const auto &v1 = vertices[vertex_indices(i, int_t(1))];
  const auto &v2 = vertices[vertex_indices(i, int_t(2))];

  if (element_type == GMSHElementType::triangle) {
    return barycenter(Triangle{v0, v1, v2});
  } else if (element_type == GMSHElementType::tetrahedron) {
    const auto &v3 = vertices[vertex_indices(i, int_t(3))];
    return barycenter(Tetrahedron{v0, v1, v2, v3});
  }

  LOG_ERR("Unknown element type.");
}

cell_centers_t compute_cell_centers(GMSHElementType element_type,
                                    const vertices_t &vertices,
                                    const vertex_indices_t &vertex_indices) {
  int_t n_cells = vertex_indices.shape(0);
  auto cell_centers = cell_centers_t(shape_t<1>{n_cells});

  for (int_t i = 0; i < n_cells; ++i) {
    cell_centers(i) = cell_center(element_type, vertices, vertex_indices, i);
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

int_t relative_neighbour_index(GMSHElementType element_type,
                               const std::vector<bool> &s) {
  return GMSHElementInfo::relative_vertex_index(element_type, s);
};

std::optional<std::pair<int_t, int_t>>
common_face(GMSHElementType element_type,
            const vertex_indices_t &vertex_indices,
            int_t i,
            int_t j) {
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

  auto ki = relative_neighbour_index(element_type, si);
  auto kj = relative_neighbour_index(element_type, sj);

  if (ki == magic_index_value) {
    return std::nullopt;
  }

  return std::optional<std::pair<int_t, int_t>>({ki, kj});
}

neighbours_t compute_neighbours(GMSHElementType element_type,
                                const vertex_indices_t &vertex_indices) {
  int_t n_vertices = count_vertices(vertex_indices);

  auto neighbours = empty_like(vertex_indices);
  std::fill(neighbours.begin(), neighbours.end(), magic_index_value);

  auto vertex_neighbours = compute_vertex_neighbours(vertex_indices);

  for (int_t vi = 0; vi < n_vertices; ++vi) {
    for (auto [i, j] : strict_symmetric_choices(vertex_neighbours[vi])) {
      auto face = common_face(element_type, vertex_indices, i, j);
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

array<Cell, 1> compute_cells(GMSHElementType element_type,
                             int_t quad_deg,
                             const vertices_t &vertices,
                             const vertex_indices_t &vertex_indices) {
  auto n_cells = vertex_indices.shape(0);
  auto cells = array<Cell, 1>(shape_t<1>{n_cells});

  if (element_type == GMSHElementType::triangle) {
    auto qr_ref = zisa::cached_triangular_quadrature_rule(quad_deg);

    for (int_t i = 0; i < n_cells; ++i) {
      const auto &v0 = vertices[vertex_indices(i, int_t(0))];
      const auto &v1 = vertices[vertex_indices(i, int_t(1))];
      const auto &v2 = vertices[vertex_indices(i, int_t(2))];

      cells[i] = Cell(zisa::denormalize(qr_ref, Triangle(v0, v1, v2)));
    }
  } else {
    auto qr_ref = zisa::make_tetrahedral_rule(quad_deg);

    for (int_t i = 0; i < n_cells; ++i) {
      const auto &v0 = vertices[vertex_indices(i, int_t(0))];
      const auto &v1 = vertices[vertex_indices(i, int_t(1))];
      const auto &v2 = vertices[vertex_indices(i, int_t(2))];
      const auto &v3 = vertices[vertex_indices(i, int_t(3))];

      cells[i] = Cell(zisa::denormalize(qr_ref, Tetrahedron(v0, v1, v2, v3)));
    }
  }

  return cells;
}

Grid::Grid(GMSHElementType element_type,
           array<XYZ, 1> vertices_,
           array<int_t, 2> vertex_indices_,
           int_t quad_deg)
    : vertex_indices(std::move(vertex_indices_)),
      vertices(std::move(vertices_)) {

  n_cells = vertex_indices.shape(0);
  n_vertices = vertices.shape(0);
  max_neighbours = vertex_indices.shape(1);

  neighbours = compute_neighbours(element_type, this->vertex_indices);
  is_valid = compute_valid_neighbours(neighbours);

  n_interior_edges = count_interior_edges(neighbours, is_valid);
  n_exterior_edges = count_exterior_edges(is_valid);
  n_edges = n_interior_edges + n_exterior_edges;

  cell_centers = compute_cell_centers(
      element_type, this->vertices, this->vertex_indices);
  edge_indices = compute_edge_indices(neighbours, is_valid);
  left_right = compute_left_right(edge_indices, neighbours, is_valid);

  volumes = compute_volumes(element_type, this->vertices, this->vertex_indices);

  normals = compute_normals(element_type,
                            this->vertices,
                            this->vertex_indices,
                            neighbours,
                            is_valid,
                            edge_indices);

  tangentials = compute_tangentials(element_type,
                                    this->normals,
                                    this->vertices,
                                    this->vertex_indices,
                                    neighbours,
                                    is_valid,
                                    edge_indices);

  if (quad_deg != 0) {
    cells = compute_cells(
        element_type, quad_deg, this->vertices, this->vertex_indices);
    normalized_moments = compute_normalized_moments(*this);
  }
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

template <class Predicate, class Ranking>
std::optional<int_t> depth_first_search(const Grid &grid,
                                        int_t i_guess,
                                        const Predicate &predicate,
                                        const Ranking &ranking) {

  if (predicate(i_guess)) {
    return i_guess;
  }

  const auto &is_valid = grid.is_valid;
  const auto &neighbours = grid.neighbours;

  auto n_cells = grid.n_cells;
  auto max_neighbours = grid.max_neighbours;

  auto visited = array<bool, 1>(shape_t<1>{n_cells});
  zisa::fill(visited, false);

  auto trail = std::stack<int_t>{};
  auto candidates = std::vector<int_t>{};
  candidates.reserve(max_neighbours);

  int_t i = i_guess;
  visited[i] = true;
  trail.push(i);

  for (int_t count = 0; count < grid.n_cells; ++count) {
    candidates.clear();

    // populate candidates.
    for (int_t k = 0; k < max_neighbours; ++k) {
      int_t cand = neighbours(i, k);
      if (is_valid(i, k) && !visited(cand)) {

        if (predicate(cand)) {
          return cand;
        }

        candidates.push_back(cand);
      }
    }

    if (!candidates.empty()) {
      auto comp
          = [&ranking](int_t i, int_t j) { return ranking(i) < ranking(j); };

      i = *std::max_element(candidates.begin(), candidates.end(), comp);

      visited[i] = true;
      trail.push(i);
    }

    else {
      if (trail.empty()) {
        return std::nullopt;
      }

      i = trail.top();
      trail.pop();
    }
  }

  return std::nullopt;
}

Triangle triangle(const Grid &grid, int_t i) {
  const auto &v0 = grid.vertex(i, int_t(0));
  const auto &v1 = grid.vertex(i, int_t(1));
  const auto &v2 = grid.vertex(i, int_t(2));

  return Triangle(v0, v1, v2);
}

Tetrahedron tetrahedron(const Grid &grid, int_t i) {
  const auto &v0 = grid.vertex(i, int_t(0));
  const auto &v1 = grid.vertex(i, int_t(1));
  const auto &v2 = grid.vertex(i, int_t(2));
  const auto &v3 = grid.vertex(i, int_t(3));

  return Tetrahedron(v0, v1, v2, v3);
}

bool is_inside_cell(const Grid &grid, int_t i, const XYZ &x) {
  GMSHElementType element_type
      = (grid.max_neighbours == 3 ? GMSHElementType::triangle
                                  : GMSHElementType::tetrahedron);

  if (element_type == GMSHElementType::triangle) {
    auto tri = triangle(grid, i);
    return is_inside(Barycentric2D(tri, x));
  } else {
    auto tet = tetrahedron(grid, i);
    return is_inside(Barycentric3D(tet, x));
  }
}

std::optional<int_t> locate(const Grid &grid, const XYZ &x, int_t i_guess) {

  auto predicate = [&grid, &x](int_t i) { return is_inside_cell(grid, i, x); };

  auto ranking = [&grid, &x](int_t i) {
    auto d = zisa::norm(grid.cell_centers(i) - x);
    return 1.0 / (1.0 + d);
  };

  return depth_first_search(grid, i_guess, predicate, ranking);
}

double volume(const Grid &grid) {
  return zisa::reduce::sum(
      zisa::cells(grid), [](int_t, const Cell &cell) { return volume(cell); });
}

std::shared_ptr<Grid> load_grid(const std::string &filename, int_t quad_deg) {

  auto len = filename.size();

  if (filename.substr(len - 4) == ".msh") {
    return load_gmsh(filename, quad_deg);
  } else if (filename.substr(len - 3) == ".h5") {
    auto reader = HDF5SerialReader(filename);
    return std::make_shared<Grid>(Grid::load(reader));
  }

  LOG_ERR(string_format("Unknown filetype. [%s]", filename.c_str()));
}

std::shared_ptr<Grid> load_gmsh(const std::string &filename, int_t quad_deg) {
  auto gmsh = GMSHData(filename);

  auto max_neighbours = GMSHElementInfo::n_vertices(gmsh.element_type);
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

  return std::make_shared<Grid>(gmsh.element_type,
                                std::move(vertices),
                                std::move(vertex_indices),
                                quad_deg);
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
  grid.tangentials = array<XYZ, 2>::load(reader, "tangentials");

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
