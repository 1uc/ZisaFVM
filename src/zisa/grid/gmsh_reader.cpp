// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <memory>

#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/utils/logging.hpp>
#include <zisa/utils/string_format.hpp>

namespace zisa {

auto GMSHElementInfo::n_vertices(GMSHElementType element_type) -> index_t {
  if (element_type == GMSHElementType::triangle) {
    return 3;
  } else if (element_type == GMSHElementType::tetrahedron) {
    return 4;
  }

  LOG_ERR(string_format("Unknown element type. [%d]", element_type));
}

GMSHElementInfo::index_t
GMSHElementInfo::relative_vertex_index(GMSHElementType element_type,
                                       GMSHElementInfo::index_t k,
                                       GMSHElementInfo::index_t rel) {
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

GMSHElementInfo::index_t
GMSHElementInfo::relative_off_vertex_index(GMSHElementType element_type,
                                           GMSHElementInfo::index_t k) {

  if (element_type == GMSHElementType::triangle) {
    return (k + 2) % 3;
  }

  if (element_type == GMSHElementType::tetrahedron) {
    if (k == 0) {
      return 2;
    } else if (k == 1) {
      return 3;
    } else if (k == 2) {
      return 1;
    } else {
      return 0;
    }
  }

  LOG_ERR("Invalid element_type.");
}

GMSHElementInfo::index_t
GMSHElementInfo::relative_vertex_index(GMSHElementType element_type,
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

} // namespace zisa
