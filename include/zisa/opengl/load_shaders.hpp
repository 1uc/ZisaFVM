#ifndef LOAD_SHADER_H_ZRCRM
#define LOAD_SHADER_H_ZRCRM

#include <string>

namespace zisa {
GLuint load_shaders(const std::string &vertex_shader,
                    const std::string &fragment_shader);

}

#endif /* end of include guard */
