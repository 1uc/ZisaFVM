// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef LOAD_SHADER_H_ZRCRM
#define LOAD_SHADER_H_ZRCRM

#if ZISA_HAS_OPENGL == 1
#include <string>

#include <zisa/config.hpp>

namespace zisa {
GLuint load_shaders(const std::string &vertex_shader,
                    const std::string &fragment_shader);

}
#endif

#endif /* end of include guard */
