#include <memory>

#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/utils/logging.hpp>
#include <zisa/utils/string_format.hpp>

namespace zisa {
/// Base class for reading GMSH files.
class GMSHReader {
protected:
  using index_t = int_t;
  using vertices_t = array<XYZ, 1>;
  using vertex_indices_t = array<int_t, 2>;

public:
  virtual ~GMSHReader() = default;
  GMSHData load(const std::string &filename) const;

protected:
  virtual GMSHElementType load_element_type() const = 0;

  virtual vertices_t load_vertices() const = 0;
  virtual vertex_indices_t
  load_elements(GMSHElementType element_type) const = 0;

protected:
  mutable std::string filename;
};

/// Reads GMSH v2 files.
class GMSHReader2 : public GMSHReader {
protected:
  virtual GMSHElementType load_element_type() const override;

  virtual vertices_t load_vertices() const override;
  virtual vertex_indices_t
  load_elements(GMSHElementType element_type) const override;

  void load_element(vertex_indices_t &vertex_indices,
                    const std::string &line,
                    GMSHElementType element_type,
                    index_t &i) const;

  index_t count_cells(GMSHElementType element_type) const;

  void
  load_vertex(vertices_t &vertices, const std::string &line, index_t i) const;
};

namespace gmsh {
std::ifstream open_msh(const std::string &filename) {
  auto is = std::ifstream(filename);
  LOG_ERR_IF(!is, string_format("Failed to open file. [%s]", filename.c_str()));

  return is;
}

std::string getline(std::ifstream &msh) {
  std::string line;
  std::getline(msh, line);
  return line;
}

bool starts_with(const std::string &line, const std::string &token) {
  return line.compare(0, token.size(), token) == 0;
}

std::ifstream &skip_to(std::ifstream &msh, const std::string &token) {
  while (msh.good() && !starts_with(getline(msh), token)) {
    continue;
  }
  return msh;
}
}

struct GMSHFileVersion {
  int major;
  double minor;
};

GMSHFileVersion gmsh_read_version(const std::string &filename) {
  auto msh = gmsh::open_msh(filename);
  gmsh::skip_to(msh, "$MeshFormat");

  auto line = std::istringstream(gmsh::getline(msh));
  double version;
  line >> version;

  return GMSHFileVersion{int(version), version - int(version)};
}

std::unique_ptr<GMSHReader> make_gmsh_reader(const GMSHFileVersion &version) {
  if (version.major == 2) {
    return std::make_unique<GMSHReader2>();
  } else if (version.major == 4) {
    LOG_ERR("GMSH file version 4.0 needs to be implemented.");
  }
  LOG_ERR("Unknown version.");
}

std::unique_ptr<GMSHReader> make_gmsh_reader(const std::string &filename) {
  auto version = gmsh_read_version(filename);
  return make_gmsh_reader(version);
}

auto GMSHElementInfo::n_vertices_per_face(GMSHElementType element_type)
    -> index_t {
  return n_vertices(element_type) - 1;
}

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

GMSHData GMSHReader::load(const std::string &filename) const {
  this->filename = filename;

  auto element_type = load_element_type();

  auto vertex_indices = load_elements(element_type);
  auto vertices = load_vertices();

  return GMSHData{std::move(vertices), std::move(vertex_indices), element_type};
}

auto GMSHReader2::load_vertices() const -> vertices_t {
  auto msh = gmsh::open_msh(filename);
  gmsh::skip_to(msh, "$Nodes");

  index_t n_vertices;
  msh >> n_vertices;
  auto discard = gmsh::getline(msh);

  auto vertices = vertices_t(n_vertices);

  for (index_t i = 0; (i < n_vertices) && msh.good(); ++i) {
    load_vertex(vertices, gmsh::getline(msh), i);
  }

  return vertices;
}

void GMSHReader2::load_vertex(vertices_t &vertices,
                              const std::string &line_,
                              index_t i) const {
  auto line = std::istringstream(line_);

  int discard;
  double x, y, z;
  line >> discard >> x >> y >> z;

  vertices(i) = XYZ{x, y, z};
}

GMSHElementType GMSHReader2::load_element_type() const {
  auto msh = gmsh::open_msh(filename);
  gmsh::skip_to(msh, "$Elements");

  index_t n_elements;
  msh >> n_elements;
  auto discard = gmsh::getline(msh);

  for (index_t i = 0; (i < n_elements) && msh.good(); ++i) {
    auto line = std::istringstream(gmsh::getline(msh));

    index_t discard, element;
    line >> discard;
    line >> element;

    if (element == index_t(GMSHElementType::tetrahedron)) {
      return GMSHElementType::tetrahedron;
    }
  }

  return GMSHElementType::triangle;
}

int_t GMSHReader2::count_cells(GMSHElementType expected_element_type) const {
  auto msh = gmsh::open_msh(filename);
  gmsh::skip_to(msh, "$Elements");

  index_t n_elements;
  msh >> n_elements;
  gmsh::getline(msh);

  int_t n_cells = 0;
  for (index_t i = 0; (i < n_elements) && msh.good(); ++i) {
    auto line = std::istringstream(gmsh::getline(msh));

    int_t discard, element_type;
    line >> discard;
    line >> element_type;
    if (element_type == index_t(expected_element_type)) {
      n_cells += 1;
    }
  }

  return n_cells;
}

auto GMSHReader2::load_elements(GMSHElementType element_type) const
    -> vertex_indices_t {

  int_t n_cells = count_cells(element_type);

  auto msh = gmsh::open_msh(filename);
  gmsh::skip_to(msh, "$Elements");

  index_t n_elements;
  msh >> n_elements;
  auto discard = gmsh::getline(msh);

  index_t max_vertices = GMSHElementInfo::n_vertices(element_type);
  auto vertex_indices = vertex_indices_t({n_cells, max_vertices});

  int_t i_cell = 0;
  for (index_t i = 0; (i < n_elements) && msh.good(); ++i) {
    load_element(vertex_indices, gmsh::getline(msh), element_type, i_cell);
  }

  return vertex_indices;
}

void GMSHReader2::load_element(vertex_indices_t &vertex_indices,
                               const std::string &line_,
                               GMSHElementType expected_element_type,
                               index_t &i) const {
  index_t n_vertices = GMSHElementInfo::n_vertices(expected_element_type);
  index_t n_tags, discard, element_type;

  auto line = std::istringstream(line_);
  line >> discard;
  line >> element_type;
  if (element_type != index_t(expected_element_type)) {
    return;
  }

  line >> n_tags;
  for (index_t k = 0; k < n_tags; ++k) {
    line >> discard;
  }

  for (index_t k = 0; k < n_vertices; ++k) {
    line >> vertex_indices(i, k);
    vertex_indices(i, k) -= 1;
  }

  i += 1;
}

std::ostream &operator<<(std::ostream &os, const GMSHData &data) {
  auto &vertex_indices = data.vertex_indices;

  os << "[ \n";
  for (int_t i = 0; i < vertex_indices.shape(0); ++i) {
    os << "  [";
    for (int_t k = 0; k < vertex_indices.shape(1); ++k) {
      os << vertex_indices(i, k);
      os << (k < vertex_indices.shape(1) - 1 ? ", " : "");
    }
    os << "]\n";
  }
  os << "]";

  return os;
}

GMSHData::GMSHData(const std::string &filename) {
  (*this) = make_gmsh_reader(filename)->load(filename);
}

GMSHData::GMSHData(array<XYZ, 1> vertices,
                   array<index_t, 2> vertex_indices,
                   GMSHElementType element_type)
    : vertices(std::move(vertices)),
      vertex_indices(std::move(vertex_indices)),
      element_type(element_type) {}

} // namespace zisa
