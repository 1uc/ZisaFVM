#include <zisa/math/cell_factory.hpp>

#include <zisa/math/cell.hpp>
#include <zisa/math/denormalized_rule.hpp>
#include <zisa/math/tetrahedral_rule.hpp>
#include <zisa/math/tetrahedron.hpp>
#include <zisa/math/triangle.hpp>
#include <zisa/math/triangular_rule.hpp>

namespace zisa {

Cell make_cell(const Triangle &tri, int_t quad_deg) {
  auto tri_qr = cached_triangular_quadrature_rule(quad_deg);
  auto qr = denormalize(tri_qr, tri);

  return Cell(std::move(qr));
}

Cell make_cell(const Tetrahedron &tet, int_t quad_deg) {
  auto tet_qr = make_tetrahedral_rule(quad_deg);
  auto qr = denormalize(tet_qr, tet);

  return Cell(std::move(qr));
}
  
}