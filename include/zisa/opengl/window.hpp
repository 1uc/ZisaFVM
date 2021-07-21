// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef WINDOW_H_KV5QN
#define WINDOW_H_KV5QN

namespace zisa {
namespace opengl {

class Window {
public:
  Window(const std::string &title, int width, int height);
  ~Window();

  GLFWwindow *ptr();

private:
  GLFWwindow *window_;
};

} // namespace opengl
} // namespace zisa
#endif /* end of include guard */
