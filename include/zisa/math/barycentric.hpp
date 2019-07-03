#ifndef BARYCENTRIC_H_3NBM5
#define BARYCENTRIC_H_3NBM5

#include <zisa/config.hpp>
#include <zisa/math/tetrahedron.hpp>
#include <zisa/math/triangle.hpp>

namespace zisa {

class Barycentric2D {
public:
  Barycentric2D() = default;
  Barycentric2D(double lambda1, double lambda2, double lambda3);
  Barycentric2D(const Triangle &tri, const XYZ &x);

  inline double operator[](int i) const { return lambda[i]; }

private:
  double lambda[3];
};

XYZ coord(const Triangle &tri, const Barycentric2D &x);
bool is_inside(const Barycentric2D &lambda);
std::ostream &operator<<(std::ostream &os, const Barycentric2D &bc);

class Barycentric3D {
public:
  Barycentric3D() = default;
  Barycentric3D(double lambda1, double lambda2, double lambda3, double lambda4);
  Barycentric3D(const Tetrahedron &tet, const XYZ &x);

  inline double operator[](int i) const { return lambda[i]; }

private:
  double lambda[4];
};

XYZ coord(const Tetrahedron &tet, const Barycentric3D &x);
bool is_inside(const Barycentric3D &lambda);
std::ostream &operator<<(std::ostream &os, const Barycentric3D &bc);

} // namespace zisa

#endif
