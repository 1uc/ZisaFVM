// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_CELL_FACTORY_HPP_CKOSE
#define ZISA_CELL_FACTORY_HPP_CKOSE

#include <zisa/config.hpp>

#include <zisa/math/cell.hpp>
#include <zisa/math/tetrahedron.hpp>
#include <zisa/math/triangle.hpp>

namespace zisa {

Cell make_cell(const Triangle &tri, int_t quad_deg);
Cell make_cell(const Tetrahedron &tet, int_t quad_deg);

}

#endif // ZISA_CELL_FACTORY_HPP
