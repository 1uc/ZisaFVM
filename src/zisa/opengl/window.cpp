#if ZISA_HAS_OPENGL == 1
#include <string>

// -- must come before all other GL stuff.
#include <GL/glew.h>
// ---------------------------------------

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include <zisa/config.hpp>
#include <zisa/opengl/window.hpp>

namespace zisa {
namespace opengl {

Window::Window(const std::string &title, int width, int height) {

  glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

  // We don't want the old OpenGL
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  // Open a window and create its OpenGL context
  window_ = glfwCreateWindow(width, height, title.c_str(), NULL, NULL);

  if (window_ == nullptr) {
    glfwTerminate();
    LOG_ERR("Failed to open a window.\n");
  }

  glfwMakeContextCurrent(window_);
  glewExperimental = true;

  if (glewInit() != GLEW_OK) {
    LOG_ERR("Failed to initialize GLEW.");
  }
}

Window::~Window() { glfwDestroyWindow(window_); }

GLFWwindow *Window::ptr() { return window_; }

} // namespace opengl
} // namespace zisa
#endif
