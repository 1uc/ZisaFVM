#ifndef GMSH_READER_H_BC5CX
#define GMSH_READER_H_BC5CX

#include <array>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <iostream>

#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

enum class GMSHElementType { triangle = 2, tetrahedron = 4 };

struct GMSHElementInfo {
private:
  using index_t = std::size_t;
  static constexpr index_t magic_index_value
      = std::numeric_limits<index_t>::max();

public:
  static index_t n_vertices(GMSHElementType element_type);
  static index_t n_vertices_per_face(GMSHElementType element_type);

  static index_t
  relative_vertex_index(GMSHElementType element_type, index_t k, index_t rel);

  static index_t relative_vertex_index(GMSHElementType element_type,
                                       std::vector<bool> s);

  static index_t relative_off_vertex_index(GMSHElementType element_type,
                                           index_t k);
};

} // namespace zisa
#endif /* end of include guard */
