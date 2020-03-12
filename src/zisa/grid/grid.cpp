#include <zisa/grid/grid.hpp>

#include <algorithm>
#include <list>
#include <map>
#include <optional>
#include <vector>

#include <zisa/config.hpp>
#include <zisa/grid/cell_range.hpp>
#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/grid/grid_decl.hpp>
#include <zisa/grid/grid_impl.hpp>
#include <zisa/grid/neighbour_range.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/io/hdf5_writer.hpp>
#include <zisa/loops/reduction/any.hpp>
#include <zisa/loops/reduction/max.hpp>
#include <zisa/loops/reduction/min.hpp>
#include <zisa/loops/reduction/sum.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/cell_factory.hpp>
#include <zisa/math/face_factory.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/math/symmetric_choices.hpp>
#include <zisa/math/tetrahedral_rule.hpp>
#include <zisa/math/tetrahedron.hpp>
#include <zisa/math/triangular_rule.hpp>
#include <zisa/utils/logging.hpp>

namespace zisa {

Triangle triangle(const vertices_t &vertices,
                  const vertex_indices_t &vertex_indices,
                  int_t i) {
  const auto &v0 = vertices(vertex_indices(i, int_t(0)));
  const auto &v1 = vertices(vertex_indices(i, int_t(1)));
  const auto &v2 = vertices(vertex_indices(i, int_t(2)));

  return Triangle(v0, v1, v2);
}

Tetrahedron tetrahedron(const vertices_t &vertices,
                        const vertex_indices_t &vertex_indices,
                        int_t i) {
  const auto &v0 = vertices(vertex_indices(i, int_t(0)));
  const auto &v1 = vertices(vertex_indices(i, int_t(1)));
  const auto &v2 = vertices(vertex_indices(i, int_t(2)));
  const auto &v3 = vertices(vertex_indices(i, int_t(3)));

  return Tetrahedron(v0, v1, v2, v3);
}

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

int_t count_edges(const neighbours_t &neighbours, const is_valid_t &is_valid) {
  return count_interior_edges(neighbours, is_valid)
         + count_exterior_edges(is_valid);
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

void enforce_standard_vertex_order(GMSHElementType element_type,
                                   vertices_t &vertices,
                                   vertex_indices_t &vertex_indices) {

  if (element_type == GMSHElementType::tetrahedron) {
    LOG_WARN("Implement this first.");
    return;
  }

  auto n_cells = vertex_indices.shape(0);
  zisa::for_each(PlainIndexRange(0, n_cells),
                 [&vertices, &vertex_indices](int_t i) {
                   auto v0 = vertices(vertex_indices(i, 0));
                   auto v1 = vertices(vertex_indices(i, 1));
                   auto v2 = vertices(vertex_indices(i, 2));

                   auto n = XYZ(zisa::cross(v1 - v0, v2 - v0));

                   if (n[2] < 0.0) {
                     std::swap(vertex_indices(i, 1), vertex_indices(i, 2));
                   }
                 });
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

  if (element_type == GMSHElementType::triangle) {
    return volume(triangle(vertices, vertex_indices, i));
  } else if (element_type == GMSHElementType::tetrahedron) {
    return volume(tetrahedron(vertices, vertex_indices, i));
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

array<double, 1> compute_inradii(GMSHElementType element_type,
                                 const vertices_t &vertices,
                                 const vertex_indices_t &vertex_indices) {
  int_t n_cells = vertex_indices.shape(0);
  auto inradii = array<double, 1>(shape_t<1>(n_cells));

  if (element_type == GMSHElementType::triangle) {
    zisa::for_each(PlainIndexRange(0, n_cells),
                   [&inradii, &vertices, &vertex_indices](int_t i) {
                     inradii(i)
                         = inradius(triangle(vertices, vertex_indices, i));
                   });
  } else if (element_type == GMSHElementType::tetrahedron) {
    zisa::for_each(PlainIndexRange(0, n_cells),
                   [&inradii, &vertices, &vertex_indices](int_t i) {
                     inradii(i)
                         = inradius(tetrahedron(vertices, vertex_indices, i));
                   });
  } else {
    LOG_ERR("Unknown element_type.");
  }

  return inradii;
}

array<double, 1> compute_circum_radii(GMSHElementType element_type,
                                      const vertices_t &vertices,
                                      const vertex_indices_t &vertex_indices) {
  int_t n_cells = vertex_indices.shape(0);
  auto circum_radii = array<double, 1>(shape_t<1>(n_cells));

  if (element_type == GMSHElementType::triangle) {
    zisa::for_each(PlainIndexRange(0, n_cells),
                   [&circum_radii, &vertices, &vertex_indices](int_t i) {
                     circum_radii(i)
                         = circum_radius(triangle(vertices, vertex_indices, i));
                   });
  } else if (element_type == GMSHElementType::tetrahedron) {
    zisa::for_each(PlainIndexRange(0, n_cells),
                   [&circum_radii, &vertices, &vertex_indices](int_t i) {
                     circum_radii(i) = circum_radius(
                         tetrahedron(vertices, vertex_indices, i));
                   });
  } else {
    LOG_ERR("Unknown element_type.");
  }

  return circum_radii;
}

XYZ cell_center(GMSHElementType element_type,
                const vertices_t &vertices,
                const vertex_indices_t &vertex_indices,
                int_t i) {
  if (element_type == GMSHElementType::triangle) {
    return barycenter(triangle(vertices, vertex_indices, i));
  } else if (element_type == GMSHElementType::tetrahedron) {
    return barycenter(tetrahedron(vertices, vertex_indices, i));
  }

  LOG_ERR("Unknown element type.");
}

template <class T>
static array<XYZ, 1> compute_barycenters(const array<T, 1> &cells) {
  auto n = cells.shape(0);
  auto centers = array<XYZ, 1>(shape_t<1>{n});

  zisa::for_each(index_range(n), [&centers, &cells](int_t i) {
    centers(i) = barycenter(cells(i));
  });

  return centers;
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
    for (int_t i = 0; i < n_cells; ++i) {
      auto tri = triangle(vertices, vertex_indices, i);
      cells[i] = make_cell(tri, quad_deg);
    }
  } else {
    for (int_t i = 0; i < n_cells; ++i) {
      auto tet = tetrahedron(vertices, vertex_indices, i);
      cells[i] = make_cell(tet, quad_deg);
    }
  }

  return cells;
}

array<Face, 1> compute_faces(GMSHElementType element_type,
                             int_t quad_deg,
                             const neighbours_t &neighbours,
                             const is_valid_t &is_valid,
                             const vertices_t &vertices,
                             const vertex_indices_t &vertex_indices,
                             const edge_indices_t &edge_indices) {

  auto n_cells = neighbours.shape(0);
  auto max_neighbours = neighbours.shape(1);
  auto n_edges = count_edges(neighbours, is_valid);

  auto faces = array<Face, 1>(shape_t<1>{n_edges});

  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < max_neighbours; ++k) {
      if (i < neighbours(i, k)) {
        if (element_type == GMSHElementType::triangle) {
          auto edge = triangle_face(vertices, vertex_indices, i, k);
          faces(edge_indices(i, k)) = make_face(edge, quad_deg);
        } else if (element_type == GMSHElementType::tetrahedron) {
          auto tri = tetrahedron_face(vertices, vertex_indices, i, k);
          faces(edge_indices(i, k)) = make_face(tri, quad_deg);
        }
      }
    }
  }

  return faces;
}

Triangle tetrahedron_face(const vertices_t &vertices,
                          const vertex_indices_t &vertex_indices,
                          int_t i,
                          int_t k) {

  auto element_type = GMSHElementType::tetrahedron;
  auto k0 = GMSHElementInfo::relative_vertex_index(element_type, k, 0);
  auto k1 = GMSHElementInfo::relative_vertex_index(element_type, k, 1);
  auto k2 = GMSHElementInfo::relative_vertex_index(element_type, k, 2);

  const auto &v0 = vertices(vertex_indices(i, k0));
  const auto &v1 = vertices(vertex_indices(i, k1));
  const auto &v2 = vertices(vertex_indices(i, k2));
  return Triangle(v0, v1, v2);
}

Edge triangle_face(const vertices_t &vertices,
                   const vertex_indices_t &vertex_indices,
                   int_t i,
                   int_t k) {
  auto element_type = GMSHElementType::triangle;
  auto k0 = GMSHElementInfo::relative_vertex_index(element_type, k, 0);
  auto k1 = GMSHElementInfo::relative_vertex_index(element_type, k, 1);

  const auto &v0 = vertices(vertex_indices(i, k0));
  const auto &v1 = vertices(vertex_indices(i, k1));

  return Edge(v0, v1);
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

  enforce_standard_vertex_order(
      element_type, this->vertices, this->vertex_indices);
  neighbours = compute_neighbours(element_type, this->vertex_indices);
  is_valid = compute_valid_neighbours(neighbours);

  n_interior_edges = count_interior_edges(neighbours, is_valid);
  n_exterior_edges = count_exterior_edges(is_valid);
  n_edges = n_interior_edges + n_exterior_edges;

  edge_indices = compute_edge_indices(neighbours, is_valid);
  left_right = compute_left_right(edge_indices, neighbours, is_valid);

  volumes = compute_volumes(element_type, this->vertices, this->vertex_indices);
  inradii = compute_inradii(element_type, this->vertices, this->vertex_indices);
  circum_radii = compute_circum_radii(
      element_type, this->vertices, this->vertex_indices);

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

    cell_centers = compute_barycenters(this->cells);

    faces = compute_faces(element_type,
                          quad_deg,
                          this->neighbours,
                          this->is_valid,
                          this->vertices,
                          this->vertex_indices,
                          this->edge_indices);

    face_centers = compute_barycenters(this->faces);

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
Face Grid::face(int_t i, int_t k) const { return faces(edge_indices(i, k)); }

double Grid::characteristic_length(int_t i) const { return circum_radii[i]; }

std::string Grid::str() const {
  double dx_min = smallest_inradius(*this);
  double dx_max = largest_circum_radius(*this);

  return string_format("n_cells : %d\n"
                       "n_vertices : %d\n"
                       "n_edges : %d\n"
                       "dx_min : %e\n"
                       "dx_max : %e\n"
                       "memory : %.2e GB\n",
                       n_cells,
                       n_vertices,
                       n_edges,
                       dx_min,
                       dx_max,
                       double(size_in_bytes()) * 1e-9);
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

std::shared_ptr<Grid> load_grid_gmsh(const std::string &filename, int_t deg) {
  deg = zisa::max(1ul, deg);

  auto gmsh = GMSHData(filename);
  return std::make_shared<Grid>(gmsh.element_type,
                                std::move(gmsh.vertices),
                                std::move(gmsh.vertex_indices),
                                deg);
}

std::shared_ptr<Grid> load_grid_gmsh_h5(const std::string &filename,
                                        int_t deg) {
  auto reader = HDF5SerialReader(filename);

  int n_dims = reader.read_scalar<int>("n_dims");
  auto element_type = (n_dims == 2 ? GMSHElementType::triangle
                                   : GMSHElementType::tetrahedron);

  auto vertex_indices = array<int_t, 2>::load(reader, "vertex_indices");
  auto vertices = array<XYZ, 1>::load(reader, "vertices");
  return std::make_shared<Grid>(
      element_type, std::move(vertices), std::move(vertex_indices), deg);
}

std::shared_ptr<Grid> load_grid(const std::string &filename, int_t quad_deg) {
  auto len = filename.size();

  if (filename.substr(len - 4) == ".msh") {
    return load_grid_gmsh(filename, quad_deg);
  } else if (filename.substr(len - 7) == ".msh.h5") {
    return load_grid_gmsh_h5(filename, quad_deg);
  } else if (filename.substr(len - 3) == ".h5") {
    auto reader = HDF5SerialReader(filename);
    return std::make_shared<Grid>(Grid::load(reader));
  }

  LOG_ERR(string_format("Unknown filetype. [%s]", filename.c_str()));
}

std::shared_ptr<Grid> load_grid(const std::string &filename) {
  return load_grid(filename, 1);
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
  save(writer, grid.inradii, "inradii");
  save(writer, grid.circum_radii, "circum_radii");
  save(writer, grid.normals, "normals");
  save(writer, grid.tangentials, "tangentials");

  writer.write_scalar(largest_circum_radius(grid), "dx_max");
  writer.write_scalar(smallest_inradius(grid), "dx_min");
}

double largest_circum_radius(const Grid &grid) {
  return zisa::reduce::max(cell_indices(grid),
                           [&grid](int_t i) { return grid.circum_radius(i); });
}

double smallest_inradius(const Grid &grid) {
  return zisa::reduce::min(cell_indices(grid),
                           [&grid](int_t i) { return grid.inradius(i); });
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
  grid.inradii = array<double, 1>::load(reader, "inradii");
  grid.circum_radii = array<double, 1>::load(reader, "circum_radii");
  grid.normals = array<XYZ, 1>::load(reader, "normals");
  grid.tangentials = array<XYZ, 2>::load(reader, "tangentials");

  grid.left_right
      = compute_left_right(grid.edge_indices, grid.neighbours, grid.is_valid);
  grid.normalized_moments = compute_normalized_moments(grid);

  return grid;
}

double Grid::inradius(int_t i) const { return inradii[i]; }
double Grid::circum_radius(int_t i) const { return circum_radii[i]; }

const XYZ &Grid::face_center(int_t i, int_t k) const {
  return face_centers(edge_indices(i, k));
}

int Grid::n_dims() const {
  // FIXME only correct if triangles and tetrahedra are the only elements.
  return (max_neighbours == 3 ? 2 : 3);
}

size_t Grid::size_in_bytes() const {
  return vertex_indices.size() * sizeof(vertex_indices[0])
         + edge_indices.size() * sizeof(edge_indices[0])
         + left_right.size() * sizeof(left_right[0])
         + neighbours.size() * sizeof(neighbours[0])
         + is_valid.size() * sizeof(is_valid[0])
         + vertices.size() * sizeof(vertices[0])
         + cell_centers.size() * sizeof(cell_centers[0])
         + face_centers.size() * sizeof(face_centers[0])
         + cells.size() * sizeof(cells[0]) + faces.size() * sizeof(faces[0])
         + volumes.size() * sizeof(volumes[0])
         + inradii.size() * sizeof(inradii[0])
         + circum_radii.size() * sizeof(circum_radii[0])
         + normals.size() * sizeof(normals[0])
         + tangentials.size() * sizeof(tangentials[0]);
}

array<double, 1>
normalized_moments(const Triangle &tri, int degree, int_t quad_deg) {
  auto length = characteristic_length(tri);
  auto cell = make_cell(tri, quad_deg);

  auto m = array<double, 1>(shape_t<1>{poly_dof<2>(degree)});
  double length_d = 1.0;
  for (int d = 0; d <= degree; ++d) {
    for (int k = 0; k <= d; ++k) {
      int l = d - k;

      m(poly_index(k, l)) = avg_moment(cell, k, l) / length_d;
    }

    length_d *= length;
  }

  return m;
}

array<double, 1>
normalized_moments(const Tetrahedron &tet, int degree, int_t quad_deg) {
  auto length = characteristic_length(tet);
  auto cell = make_cell(tet, quad_deg);

  auto moments = array<double, 1>(shape_t<1>{poly_dof<2>(degree)});
  double length_d = 1.0;
  for (int d = 0; d <= degree; ++d) {
    for (int k = 0; k <= d; ++k) {
      for (int l = 0; l <= d - k; ++l) {
        int m = d - k - l;

        moments(poly_index(k, l, m)) = avg_moment(cell, k, l, m) / length_d;
      }
    }

    length_d *= length;
  }

  return moments;
}

bool is_boundary_cell(const Grid &grid, int_t i) {
  auto is_boundary_edge = [&grid, i](int_t k) { return !grid.is_valid(i, k); };

  return zisa::reduce::any(
      serial_policy{}, zisa::neighbour_index_range(grid), is_boundary_edge);
}

int_t distance_to_boundary(const Grid &grid, int_t i, int_t max_distance) {
  if (max_distance == 1) {
    return is_boundary_cell(grid, i) ? 0 : max_distance;
  }

  auto distance = [&grid, i, max_distance](int_t k) {
    return 1
           + distance_to_boundary(
               grid, grid.neighbours(i, k), max_distance - 1);
  };

  return zisa::reduce::min(
      serial_policy{}, zisa::neighbour_index_range(grid), distance);
}

} // namespace zisa
