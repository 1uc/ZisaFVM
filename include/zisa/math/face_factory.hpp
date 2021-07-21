// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_FACE_FACTORY_HPP_KEOWV
#define ZISA_FACE_FACTORY_HPP_KEOWV

#include <zisa/config.hpp>

#include <zisa/math/edge.hpp>
#include <zisa/math/face.hpp>
#include <zisa/math/triangle.hpp>

namespace zisa {

Face make_face(const Edge &edge, int_t quad_deg);
Face make_face(const Triangle &tri, int_t quad_deg);

}

#endif // ZISA_FACE_FACTORY_HPP
