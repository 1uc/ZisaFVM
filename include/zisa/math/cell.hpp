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

inline bool operator==(const Cell &a, const Cell &b) { return a.qr == b.qr; }
inline bool operator!=(const Cell &a, const Cell &b) { return !(a == b); }

inline double volume(const Cell &cell) { return volume(cell.qr); }

}

#endif /* end of include guard */
