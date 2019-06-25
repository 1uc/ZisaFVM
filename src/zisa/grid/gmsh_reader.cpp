#include <memory>

#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/utils/logging.hpp>
#include <zisa/utils/string_format.hpp>

namespace zisa {
class GMSHReader {
public:
  virtual ~GMSHReader() = default;
  virtual GMSHData load(const std::string &filename) const = 0;

protected:
  mutable std::string filename;
};

class GMSHReader2 : public GMSHReader {
private:
  using index_t = std::size_t;
  using vertices_t = std::vector<std::array<double, 3>>;
  using vertex_indices_t = std::vector<std::vector<index_t>>;

public:
  virtual GMSHData load(const std::string &filename) const override;

private:
  GMSHElementType load_element_type() const;

  vertex_indices_t load_elements(GMSHElementType element_type) const;
  void load_element(vertex_indices_t &vertex_indices,
                    const std::string &line,
                    GMSHElementType element_type) const;

  void load_vertex(vertices_t &vertices, const std::string &line) const;
  vertices_t load_vertices() const;
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

std::unique_ptr<GMSHReader> gmsh_reader(const GMSHFileVersion &version) {
  if (version.major == 2) {
    return std::make_unique<GMSHReader2>();
  } else if (version.major == 4) {
    LOG_ERR("Needs to be implemented.");
  }
  LOG_ERR("Unknown version.");
}

std::unique_ptr<GMSHReader> gmsh_reader(const std::string &filename) {
  auto version = gmsh_read_version(filename);
  return gmsh_reader(version);
}

GMSHData GMSHReader2::load(const std::string &filename) const {
  this->filename = filename;

  auto element_type = load_element_type();

  auto vertex_indices = load_elements(element_type);
  auto vertices = load_vertices();

  return GMSHData{vertices, vertex_indices, element_type};
}

auto GMSHReader2::load_vertices() const -> vertices_t {
  auto vertices = vertices_t{};

  auto msh = gmsh::open_msh(filename);
  gmsh::skip_to(msh, "$Nodes");

  index_t n_vertices;
  msh >> n_vertices;
  auto discard = gmsh::getline(msh);

  vertices.reserve(n_vertices);

  for (index_t i = 0; (i < n_vertices) && msh.good(); ++i) {
    load_vertex(vertices, gmsh::getline(msh));
  }

  return vertices;
}

void GMSHReader2::load_vertex(vertices_t &vertices,
                              const std::string &line_) const {
  auto line = std::istringstream(line_);

  int discard;
  double x, y, z;
  line >> discard >> x >> y >> z;

  vertices.push_back({x, y, z});
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

auto GMSHReader2::load_elements(GMSHElementType element_kind) const
    -> vertex_indices_t {
  auto vertex_indices = vertex_indices_t{};

  auto msh = gmsh::open_msh(filename);
  gmsh::skip_to(msh, "$Elements");

  index_t n_elements;
  msh >> n_elements;
  auto discard = gmsh::getline(msh);

  vertex_indices.reserve(n_elements);

  for (index_t i = 0; (i < n_elements) && msh.good(); ++i) {
    load_element(vertex_indices, gmsh::getline(msh), element_kind);
  }

  return vertex_indices;
}

auto GMSHElementInfo::n_vertices(GMSHElementType element_type) -> index_t {
  if (element_type == GMSHElementType::triangle) {
    return 3;
  } else if (element_type == GMSHElementType::tetrahedron) {
    return 4;
  }

  LOG_ERR(string_format("Unknown element type. [%d]", element_type));
}

void GMSHReader2::load_element(vertex_indices_t &vertex_indices,
                               const std::string &line_,
                               GMSHElementType element_type) const {
  index_t n_vertices = GMSHElementInfo::n_vertices(element_type);
  index_t n_tags, discard, element;

  auto line = std::istringstream(line_);
  line >> discard;
  line >> element;
  if (element != index_t(element_type)) {
    return;
  }

  line >> n_tags;
  for (index_t i = 0; i < n_tags; ++i) {
    line >> discard;
  }

  vertex_indices.emplace_back(n_vertices);
  auto &triangle = vertex_indices.back();

  for (index_t i = 0; i < n_vertices; ++i) {
    line >> triangle[i];
    --triangle[i];
  }
}

std::ostream &operator<<(std::ostream &os, const GMSHData &data) {
  auto &vertex_indices = data.vertex_indices;

  os << "[ \n";
  for (auto &&triangle : vertex_indices) {
    os << "  [" << triangle[0] << ", " << triangle[1] << ", " << triangle[2]
       << "]\n";
  }
  os << "]";

  return os;
}

GMSHData::GMSHData(const std::string &filename) {
  (*this) = gmsh_reader(filename)->load(filename);
}

GMSHData::GMSHData(std::vector<std::array<double, 3>> vertices,
                   std::vector<std::vector<GMSHData::index_t>> vertex_indices,
                   GMSHElementType element_type)
    : vertices(std::move(vertices)),
      vertex_indices(std::move(vertex_indices)),
      element_type(element_type) {}

} // namespace zisa
