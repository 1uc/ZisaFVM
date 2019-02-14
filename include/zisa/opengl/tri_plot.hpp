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

#include <zisa/io/colors.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/opengl/window.hpp>

namespace zisa {
namespace opengl {

class TriPlot {
private:
  using Vertices = zisa::array<XYZ, 1>;
  using VertexIndices = zisa::array<int_t, 2>;

public:
  TriPlot(std::shared_ptr<Window> window,
          const array<XYZ, 1> &vertices,
          const array<int_t, 2> &vertex_indices);

  void draw(const array<RGBColor, 1> &colors) const;

private:
  void init_vao();
  void init_vbos();
  void init_shaders();

  void compute_normalized_coordinates(const Vertices &vertices,
                                      const VertexIndices &vertex_indices);

  void clear() const;

  void bind_vertices() const;
  void bind_colors(const array<RGBColor, 1> &colors) const;

private:
  std::vector<std::array<float, 3>> vertices;

  // -- Internal OpenGL stuff
  std::shared_ptr<Window> window;
  GLuint program_id;

  GLuint vertex_buffer;
  GLuint color_buffer;
};

} // namespace opengl
} // namespace zisa
#endif

#endif /* end of include guard */
