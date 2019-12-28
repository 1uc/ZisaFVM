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
};

struct GMSHData {
private:
  using index_t = std::size_t;

public:
  std::vector<std::array<double, 3>> vertices;
  std::vector<std::vector<index_t>> vertex_indices;
  GMSHElementType element_type;

public:
  GMSHData(const std::string &filename);
  GMSHData(std::vector<std::array<double, 3>> vertices,
           std::vector<std::vector<index_t>> vertex_indices,
           GMSHElementType element_type);
};

std::ostream &operator<<(std::ostream &os, const GMSHData &gmsh_data);

} // namespace zisa
#endif /* end of include guard */
