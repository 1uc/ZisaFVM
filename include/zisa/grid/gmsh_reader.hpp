#ifndef GMSH_READER_H_BC5CX
#define GMSH_READER_H_BC5CX

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <iostream>

namespace zisa {

struct GMSHElementInfo {
private:
  using index_t = std::size_t;

public:
  static index_t n_vertices(index_t element_type);
};

struct GMSHData {
private:
  using index_t = std::size_t;

public:
  std::vector<std::array<double, 3>> vertices;
  std::vector<std::vector<index_t>> vertex_indices;

public:
  GMSHData(const std::string &filename);
  GMSHData(std::vector<std::array<double, 3>> vertices,
           std::vector<std::vector<index_t>> vertex_indices);
};

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
  //  void load_triangles();
  //  void load_tetrahedra();

  vertex_indices_t load_elements(index_t element_kind) const;
  void load_element(vertex_indices_t &vertex_indices,
                    const std::string &line,
                    index_t element_kind) const;

  void load_vertex(vertices_t &vertices, const std::string &line) const;
  vertices_t load_vertices() const;
};

std::ostream &operator<<(std::ostream &os, const GMSHReader &reader);

} // namespace zisa
#endif /* end of include guard */
