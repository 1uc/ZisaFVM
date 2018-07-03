#ifndef GMSH_READER_H_BC5CX
#define GMSH_READER_H_BC5CX

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <iostream>

namespace zisa {

class GMSHReader {
private:
  using index_t = std::size_t;

public:
  std::vector<std::vector<double>> vertices;
  std::vector<std::vector<index_t>> vertex_indices;

public:
  GMSHReader(const std::string &filename);

private:
  std::ifstream open_msh() const;

  void load_triangles();
  void load_elements(index_t element_kind);
  void load_element(std::istringstream line, index_t element_kind);
  void load_vertex(std::istringstream line);
  void load_vertices();

  std::string getline(std::ifstream &msh) const;
  std::ifstream &skip_to(std::ifstream &msh, const std::string &token) const;
  bool starts_with(const std::string &line, const std::string &token) const;

private:
  std::string filename;
};

std::ostream &operator<<(std::ostream &os, const GMSHReader &reader);

} // namespace zisa
#endif /* end of include guard */
