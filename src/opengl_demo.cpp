#include <random>

#include <zisa/grid/grid.hpp>
#include <zisa/opengl/tri_plot.hpp>

int main() {
#if ZISA_HAS_OPENGL == 1
  glewExperimental = true;
  if (!glfwInit()) {
    std::cerr << "Failed to initialize GLFW.\n";
    return -1;
  }

  std::random_device rd;
  std::mt19937 e(rd());
  std::uniform_real_distribution<float> dist(0, 1.0);

  auto window = std::make_shared<zisa::Window>("OpenGL demo", 1024, 1024);

  auto grid = zisa::load_gmsh("grids/convergence/unit_square_5.msh");
  zisa::TriPlot plot(window, grid->vertices, grid->vertex_indices);

  auto colors = std::vector<std::array<float, 3>>(grid->n_cells * 3);

  // Ensure we can capture the escape key being pressed below
  glfwSetInputMode(window->ptr(), GLFW_STICKY_KEYS, GL_TRUE);

  do {
    for (zisa::int_t i = 0; i < colors.size() / 3; ++i) {
      colors[3*i][0] = dist(e);
      colors[3*i][1] = dist(e);
      colors[3*i][2] = dist(e);

      colors[3*i + 1][0] = colors[3*i][0];
      colors[3*i + 1][1] = colors[3*i][1];
      colors[3*i + 1][2] = colors[3*i][2];

      colors[3*i + 2][0] = colors[3*i][0];
      colors[3*i + 2][1] = colors[3*i][1];
      colors[3*i + 2][2] = colors[3*i][2];
    }
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
