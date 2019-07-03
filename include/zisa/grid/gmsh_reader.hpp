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

  static index_t
  relative_vertex_index(GMSHElementType element_type, index_t k, index_t rel) {
    if (element_type == GMSHElementType::triangle) {
      return (k + rel) % 3;
    }

    if (element_type == GMSHElementType::tetrahedron) {
      // clang-format off
      if (k == 0) {
        if (rel == 0)      { return 0; }
        else if (rel == 1) { return 1; }
        else               { return 3; }
      }
      else if (k == 1) {
        if (rel == 0)      { return 0; }
        else if (rel == 1) { return 2; }
        else               { return 1; }
      }
      else if (k == 2) {
        if (rel == 0)      { return 0; }
        else if (rel == 1) { return 3; }
        else               { return 2; }
      }
      else {
        if (rel == 0)      { return 1; }
        else if (rel == 1) { return 2; }
        else               { return 3; }
      }
      // clang-format on
    }

    LOG_ERR("Unknown element type.");
  }

  static index_t relative_vertex_index(GMSHElementType element_type,
                                       std::vector<bool> s) {

    if (element_type == GMSHElementType::triangle) {
      for (index_t k = 0; k < 3; ++k) {
        if (s[k] && s[(k + 1) % 3]) {
          return k;
        }
      }
      return magic_index_value;
    } else if (element_type == GMSHElementType::tetrahedron) {

      if (s[0] && s[1] && s[3]) {
        return index_t(0);
      } else if (s[0] && s[2] && s[1]) {
        return index_t(1);
      } else if (s[0] && s[3] && s[2]) {
        return index_t(2);
      } else if (s[1] && s[2] && s[3]) {
        return index_t(3);
      } else {
        return magic_index_value;
      }
    }
    return magic_index_value;
  }
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
