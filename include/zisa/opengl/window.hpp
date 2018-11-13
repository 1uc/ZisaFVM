#ifndef WINDOW_H_KV5QN
#define WINDOW_H_KV5QN

namespace zisa {

class Window {
public:
  Window(const std::string &title, int width, int height);
  ~Window();

  GLFWwindow *ptr();

private:
  GLFWwindow *window_;
};

} // namespace zisa
#endif /* end of include guard */
