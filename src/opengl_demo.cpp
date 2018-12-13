#include <random>

#include <zisa/grid/grid.hpp>
#include <zisa/opengl/tri_plot.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/io/colors.hpp>
#include <zisa/io/color_map.hpp>

int main() {
#if ZISA_HAS_OPENGL == 1
  glewExperimental = true;
  if (!glfwInit()) {
    std::cerr << "Failed to initialize GLFW.\n";
    return -1;
  }

  auto window = std::make_shared<zisa::opengl::Window>("OpenGL demo", 600, 600);

  auto grid = zisa::load_gmsh("grids/convergence/unit_square_2.msh");
  zisa::opengl::TriPlot plot(window, grid->vertices, grid->vertex_indices);

  auto colors = zisa::array<zisa::RGBColor, 1>(zisa::shape_t<1>{grid->n_cells * 3});

  // Ensure we can capture the escape key being pressed below
  glfwSetInputMode(window->ptr(), GLFW_STICKY_KEYS, GL_TRUE);

  // auto blue = zisa::LABColor{35.7f, -8.4f, -27.2f};
  auto blue_ = zisa::RGBColor{0.0f, 0.353f, 0.498f};
  // auto blue__ = zisa::XYZColor{7.5f, 8.8f, 21.4f};

  // INFO(xyz2rgb(blue__));
  // REQUIRE(almost_equal(xyz2rgb(blue__), blue_, 0.1));

  // auto red = zisa::LABColor{30.4f, 37.2f, 15.7f};
  // auto red_ = zisa::RGBColor{0.498f, 0.165f, 0.192f};

  auto colormap = zisa::ColorMap(11);

  for (zisa::int_t i = 0; i < colors.size() / 3; ++i) {
    auto n_cells = colors.size() / 3;

    colors[3 * i] = colormap(double(i) / n_cells * 2.0 - 1.0);

    colors[3 * i + 1] = colors[3 * i];
    colors[3 * i + 2] = colors[3 * i];
  }

  do {
    std::cout << "." << std::flush;
    plot.draw(colors);

    glfwPollEvents();

  } while (glfwGetKey(window->ptr(), GLFW_KEY_ESCAPE) != GLFW_PRESS
           && !glfwWindowShouldClose(window->ptr()));

  return 0;

#else
  std::cerr << "Compiled without OpenGL support. \n";
  return -1;
#endif
}
