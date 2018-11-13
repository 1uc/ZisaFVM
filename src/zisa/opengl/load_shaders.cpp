#if ZISA_HAS_OPENGL == 1

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

// -- must come before all other GL stuff.
#include <GL/glew.h>
// ---------------------------------------

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include <zisa/opengl/load_shaders.hpp>

namespace zisa {
static std::string read_file(char const *const file_name) {

  std::ifstream file(file_name, std::ios::in);
  if (file.good()) {
    std::stringstream sstr;
    sstr << file.rdbuf();
    return sstr.str();
  }

  std::cerr << "Failed to load '" << file_name << "'\n";
  std::exit(-1);

  return "";
}

static void check_gl_shader_error(GLuint id) {
  int log_length;
  glGetShaderiv(id, GL_INFO_LOG_LENGTH, &log_length);
  if (log_length > 0) {
    std::vector<char> error_msg(static_cast<int_t>(log_length + 1));
    glGetShaderInfoLog(id, log_length, nullptr, error_msg.data());

    printf("%s\n", error_msg.data());
    std::exit(-1);
  }
}

static void check_gl_program_error(GLuint id) {
  int log_length;
  glGetProgramiv(id, GL_INFO_LOG_LENGTH, &log_length);
  if (log_length > 0) {
    std::vector<char> error_msg(static_cast<int_t>(log_length + 1));
    glGetProgramInfoLog(id, log_length, nullptr, error_msg.data());

    printf("%s\n", error_msg.data());
    std::exit(-1);
  }
}

static GLuint compile_shader(const std::string &code, GLuint shader_type) {
  GLuint id = glCreateShader(shader_type);

  char const *code_ptr = code.c_str();
  glShaderSource(id, 1, &code_ptr, nullptr);
  glCompileShader(id);

  check_gl_shader_error(id);

  return id;
}

static GLuint link_shaders(const std::vector<GLuint> &shaders) {
  GLuint program_id = glCreateProgram();

  for (auto &s : shaders) {
    glAttachShader(program_id, s);
  }
  glLinkProgram(program_id);
  check_gl_program_error(program_id);

  for (auto &s : shaders) {
    glDetachShader(program_id, s);
  }

  return program_id;
}

static GLuint load_shaders(char const *const vertex_file_path,
                           char const *const fragment_file_path) {

  auto vertex_code = read_file(vertex_file_path);
  auto fragment_code = read_file(fragment_file_path);

  auto vertex_id = compile_shader(vertex_code, GL_VERTEX_SHADER);
  auto fragment_id = compile_shader(fragment_code, GL_FRAGMENT_SHADER);

  auto program_id = link_shaders({vertex_id, fragment_id});

  glDeleteShader(vertex_id);
  glDeleteShader(fragment_id);

  return program_id;
}

GLuint load_shaders(const std::string &vertex_file_path,
                    const std::string &fragment_file_path) {

  return load_shaders(vertex_file_path.c_str(), fragment_file_path.c_str());
}

} // namespace zisa

#endif
