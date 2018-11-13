#ifndef TRI_PLOT_H_Q21GL
#define TRI_PLOT_H_Q21GL

#if ZISA_HAS_OPENGL == 1

#include <array>
#include <memory>
#include <string>
#include <vector>

// -- must come before all other GL stuff.
#include <GL/glew.h>
// ---------------------------------------

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include <zisa/math/cartesian.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/opengl/window.hpp>

namespace zisa {

class TriPlot {
private:
  using Vertices = zisa::array<XY, 1>;
  using VertexIndices = zisa::array<int_t, 2>;

public:
  TriPlot(std::shared_ptr<Window> window,
          const Vertices &vertices,
          const VertexIndices &vertex_indices);

  void draw(const std::vector<std::array<float, 3>> &colors) const;

private:
  void init_vao();
  void init_vbos();
  void init_shaders();

  void compute_normalized_coordinates(const Vertices &vertices,
                                      const VertexIndices &vertex_indices);

  void clear() const;

  void bind_vertices() const;
  void bind_colors(const std::vector<std::array<float, 3>> &colors) const;

private:
  std::vector<std::array<float, 3>> vertices;

  // -- Internal OpenGL stuff
  std::shared_ptr<Window> window;
  GLuint program_id;

  GLuint vertex_buffer;
  GLuint color_buffer;
};

} // namespace zisa
#endif

#endif /* end of include guard */
