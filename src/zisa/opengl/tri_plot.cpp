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

#include <zisa/opengl/load_shaders.hpp>
#include <zisa/opengl/tri_plot.hpp>
#include <zisa/opengl/window.hpp>

namespace zisa {

TriPlot::TriPlot(std::shared_ptr<Window> window,
                 const Vertices &vertices,
                 const VertexIndices &vertex_indices)

    : vertices(vertex_indices.size()), window(window) {

  assert(window != nullptr);

  compute_normalized_coordinates(vertices, vertex_indices);

  init_vao();
  init_vbos();
  init_shaders();
}

void TriPlot::init_vao() {
  GLuint vertex_array;
  glGenVertexArrays(1, &vertex_array);
  glBindVertexArray(vertex_array);
}

void TriPlot::init_vbos() {
  std::size_t vertex_buffer_size = vertices.size() * sizeof(vertices[0]);

  glGenBuffers(1, &vertex_buffer);
  glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
  glBufferData(GL_ARRAY_BUFFER,
               vertex_buffer_size,
               (void *)&vertices[0],
               GL_STATIC_DRAW);

  glGenBuffers(1, &color_buffer);
}

void TriPlot::init_shaders() {
  program_id = load_shaders("shaders/vertex_shader.glsl",
                            "shaders/fragment_shader.glsl");
  glUseProgram(program_id);
}

void TriPlot::compute_normalized_coordinates(
    const Vertices &vertices, const VertexIndices &vertex_indices) {

  std::array<float, 2> min_, max_;
  for (int_t k = 0; k < 2; ++k) {
    min_[k] = std::numeric_limits<float>::max();
    max_[k] = std::numeric_limits<float>::min();
  }

  for (const auto &v : vertices) {
    for (int_t k = 0; k < 2; ++k) {
      min_[k] = std::min(min_[k], float(v[k]));
      max_[k] = std::max(max_[k], float(v[k]));
    }
  }

  auto n_cells = vertex_indices.shape(0);
  for (int_t i = 0; i < n_cells; ++i) {
    for (int_t k = 0; k < 3; ++k) {
      int_t j = vertex_indices(i, k); // index in 'vertices'.
      int_t jj = 3 * i + k;           // index in 'this->vertices'.

      for (int_t kk = 0; kk < 2; ++kk) {
        this->vertices[jj][kk]
            = float((vertices[j][kk] - min_[kk]) / (max_[kk] - min_[kk]));
      }
      this->vertices[jj][2] = 0.0;
    }
  }
}

void TriPlot::clear() const {
  glClearColor(1, 1, 1, 1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void TriPlot::bind_vertices() const {
  glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void *)0);
}

void TriPlot::bind_colors(
    const std::vector<std::array<float, 3>> &colors) const {
  glBindBuffer(GL_ARRAY_BUFFER, color_buffer);
  glBufferData(GL_ARRAY_BUFFER,
               colors.size() * sizeof(colors[0]),
               (void *)colors.data(),
               GL_STREAM_DRAW);

  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void *)0);
}

void TriPlot::draw(const std::vector<std::array<float, 3>> &colors) const {
  assert(colors.size() == vertices.size());

  clear();

  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);

  bind_vertices();
  bind_colors(colors);

  glDrawArrays(GL_TRIANGLES, 0, int(colors.size()));

  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);

  glfwSwapBuffers(window->ptr());
}

} // namespace zisa

#endif
