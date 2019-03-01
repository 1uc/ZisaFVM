#ifndef BARYCENTRIC_H_3NBM5
#define BARYCENTRIC_H_3NBM5

#include <zisa/config.hpp>
#include <zisa/math/triangle.hpp>

namespace zisa {

class Barycentric {
public:
  Barycentric() = default;
  Barycentric(double lambda1, double lambda2, double lambda3);
  Barycentric(const Triangle &tri, const XYZ &x);

  double operator[](int i) const { return lambda[i]; }

private:
  double lambda[3];
};

XYZ coord(const Triangle &tri, const Barycentric &x);

std::ostream &operator<<(std::ostream &os, const Barycentric &bc);

} // namespace zisa

#endif
