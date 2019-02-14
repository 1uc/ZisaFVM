#include <zisa/io/tri_plot.hpp>

#if ZISA_HAS_OPENGL != 0
#include <zisa/opengl/window.hpp>
#endif

namespace zisa {
TriPlot::TriPlot(const array<XYZ, 1> &vertices,
                 const array<int_t, 2> &vertex_indices,
                 const std::string &title) {

#if ZISA_HAS_OPENGL == 0
  LOG_ERR("Compiled without OpenGL support.");
  ZISA_UNUSED(vertices);
  ZISA_UNUSED(vertex_indices);
  ZISA_UNUSED(title);
#else

  auto window = std::make_shared<opengl::Window>(title, 600, 600);
  tri_plot_ = std::make_unique<opengl::TriPlot>(std::move(window), vertices, vertex_indices);

#endif
}

void TriPlot::draw(const array<RGBColor, 1> &colors) const {
#if ZISA_HAS_OPENGL == 0
  LOG_ERR("Compiled without OpenGL support.");
  ZISA_UNUSED(colors);
#else
  tri_plot_->draw(colors);
#endif
}

} // namespace zisa
