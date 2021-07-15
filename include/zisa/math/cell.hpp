// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef CELL_H_FC29I
#define CELL_H_FC29I

#include <zisa/config.hpp>

#include <zisa/math/denormalized_rule.hpp>

namespace zisa {
class Cell {
public:
  DenormalizedRule qr;

  Cell() = default;
  explicit Cell(DenormalizedRule qr) : qr(std::move(qr)) {}
};

std::string str(const Cell &cell);

inline bool operator==(const Cell &a, const Cell &b) { return a.qr == b.qr; }
inline bool operator!=(const Cell &a, const Cell &b) { return !(a == b); }

XYZ barycenter(const Cell &cell);
inline double volume(const Cell &cell) { return volume(cell.qr); }

double avg_moment(const Cell &cell, int x_deg, int y_deg);
double avg_moment(const Cell &cell, int x_deg, int y_deg, int z_deg);

}

#endif /* end of include guard */
