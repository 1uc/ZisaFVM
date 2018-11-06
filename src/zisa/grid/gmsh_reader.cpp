#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/utils/logging.hpp>
#include <zisa/utils/string_format.hpp>

namespace zisa {

GMSHReader::GMSHReader(const std::string &filename) : filename(filename) {
  load_vertices();
  load_triangles();
}

void GMSHReader::load_vertices() {
  auto msh = open_msh();
  skip_to(msh, "$Nodes");

  index_t n_vertices;
  msh >> n_vertices;
  auto discard = getline(msh);

  vertices.reserve(n_vertices);

  for (index_t i = 0; (i < n_vertices) && msh.good(); ++i) {
    load_vertex(std::istringstream(getline(msh)));
  }
}

void GMSHReader::load_vertex(std::istringstream line) {
  int discard;
  double x, y, z;
  line >> discard >> x >> y >> z;

  vertices.push_back({x, y, z});
}

void GMSHReader::load_triangles() { load_elements(2); }

void GMSHReader::load_elements(index_t element_kind) {

  auto msh = open_msh();
  skip_to(msh, "$Elements");

  index_t n_elements;
  msh >> n_elements;
  auto discard = getline(msh);

  vertex_indices.reserve(n_elements);

  for (index_t i = 0; (i < n_elements) && msh.good(); ++i) {
    load_element(std::istringstream(getline(msh)), element_kind);
  }
}

void GMSHReader::load_element(std::istringstream line, index_t element_kind) {
  assert(element_kind == 2);

  index_t n_vertices = 3;
  index_t n_tags, discard, element;

  line >> discard;
  line >> element;
  if (element != element_kind) {
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

std::string GMSHReader::getline(std::ifstream &msh) const {
  std::string line;
  std::getline(msh, line);
  return line;
}

std::ifstream GMSHReader::open_msh() const {
  auto is = std::ifstream(filename);
  LOG_ERR_IF(!is, string_format("Failed to open file. [%s]", filename.c_str()));

  return is;
}

std::ifstream &GMSHReader::skip_to(std::ifstream &msh,
                                   const std::string &token) const {
  while (msh.good() && !starts_with(getline(msh), token)) {
    continue;
  }
  return msh;
}

bool GMSHReader::starts_with(const std::string &line,
                             const std::string &token) const {
  return line.compare(0, token.size(), token) == 0;
}

std::ostream &operator<<(std::ostream &os, const GMSHReader &reader) {
  auto &vertex_indices = reader.vertex_indices;

  os << "[ \n";
  for (auto &&triangle : vertex_indices) {
    os << "  [" << triangle[0] << ", " << triangle[1] << ", " << triangle[2]
       << "]\n";
  }
  os << "]";

  return os;
}

} // namespace zisa
