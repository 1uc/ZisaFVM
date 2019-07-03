/* Detail functions to compute the grid.
 */

#ifndef GRID_IMPL_H_8WPAM
#define GRID_IMPL_H_8WPAM

#include <limits>
#include <map>

#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/grid/grid_decl.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {
constexpr int_t magic_index_value = std::numeric_limits<int_t>::max();

using vertices_t = array<XYZ, 1>;
using vertex_indices_t = array<int_t, 2>;
using vertex_neighbours_t = array<std::vector<int_t>, 1>;
using neighbours_t = array<int_t, 2>;
using edge_indices_t = array<int_t, 2>;
using left_right_t = array<std::pair<int_t, int_t>, 1>;
using volumes_t = array<double, 1>;
using normals_t = array<XYZ, 1>;
using tangentials_t = array<XYZ, 2>;
using cell_centers_t = array<XYZ, 1>;
using is_valid_t = array<bool, 2>;

neighbours_t compute_neighbours(GMSHElementType element_type,
                                const vertex_indices_t &vertex_indices);
is_valid_t compute_valid_neighbours(const neighbours_t &neighbours);
edge_indices_t compute_edge_indices(const neighbours_t &neighbours,
                                    const is_valid_t &is_valid);
left_right_t compute_left_right(const edge_indices_t &edge_indices,
                                const neighbours_t &neighbours,
                                const is_valid_t &is_valid);

normals_t compute_normals(GMSHElementType element_type,
                          const vertices_t &vertices,
                          const vertex_indices_t &vertex_indices,
                          const neighbours_t &neighbours,
                          const is_valid_t &is_valid,
                          const edge_indices_t &edge_indices);

volumes_t compute_volumes(GMSHElementType element_type,
                          const vertices_t &vertices,
                          const vertex_indices_t &vertex_indices);
} // namespace zisa
#endif /* end of include guard */
