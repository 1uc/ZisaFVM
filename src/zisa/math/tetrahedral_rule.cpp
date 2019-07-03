#include <zisa/math/tetrahedral_rule.hpp>

namespace zisa {

TetrahedralRule make_tetrahedral_rule(int_t deg) {
  if (deg <= 1) {
    auto qr = TetrahedralRule(1);

    qr.weights[0] = 1.0;
    qr.points[0] = Barycentric3D{0.25, 0.25, 0.25, 0.25};

    return qr;
  }

  if (deg <= 3) {
    auto qr = TetrahedralRule(4);

    qr.weights[0] = 0.25;
    qr.points[0] = Barycentric3D{0.5854101966249680,
                                 0.1381966011250110,
                                 0.1381966011250110,
                                 0.1381966011250110};

    qr.weights[1] = 0.25;
    qr.points[1] = Barycentric3D{0.1381966011250110,
                                 0.5854101966249680,
                                 0.1381966011250110,
                                 0.1381966011250110};

    qr.weights[2] = 0.25;
    qr.points[2] = Barycentric3D{0.1381966011250110,
                                 0.1381966011250110,
                                 0.5854101966249680,
                                 0.1381966011250110};

    qr.weights[3] = 0.25;
    qr.points[3] = Barycentric3D{0.1381966011250110,
                                 0.1381966011250110,
                                 0.1381966011250110,
                                 0.5854101966249680};

    return qr;
  }

  if (deg <= 5) {
    auto qr = TetrahedralRule(10);

    qr.weights[0] = 0.0476331348432089;
    qr.points[0] = Barycentric3D{0.7784952948213300,
                                 0.0738349017262234,
                                 0.0738349017262234,
                                 0.0738349017262234};

    qr.weights[1] = 0.0476331348432089;
    qr.points[1] = Barycentric3D{0.0738349017262234,
                                 0.7784952948213300,
                                 0.0738349017262234,
                                 0.0738349017262234};

    qr.weights[2] = 0.0476331348432089;
    qr.points[2] = Barycentric3D{0.0738349017262234,
                                 0.0738349017262234,
                                 0.7784952948213300,
                                 0.0738349017262234};

    qr.weights[3] = 0.0476331348432089;
    qr.points[3] = Barycentric3D{0.0738349017262234,
                                 0.0738349017262234,
                                 0.0738349017262234,
                                 0.7784952948213300};

    qr.weights[4] = 0.1349112434378610;
    qr.points[4] = Barycentric3D{0.4062443438840510,
                                 0.4062443438840510,
                                 0.0937556561159491,
                                 0.0937556561159491};

    qr.weights[5] = 0.1349112434378610;
    qr.points[5] = Barycentric3D{0.4062443438840510,
                                 0.0937556561159491,
                                 0.4062443438840510,
                                 0.0937556561159491};

    qr.weights[6] = 0.1349112434378610;
    qr.points[6] = Barycentric3D{0.4062443438840510,
                                 0.0937556561159491,
                                 0.0937556561159491,
                                 0.4062443438840510};

    qr.weights[7] = 0.1349112434378610;
    qr.points[7] = Barycentric3D{0.0937556561159491,
                                 0.4062443438840510,
                                 0.4062443438840510,
                                 0.0937556561159491};

    qr.weights[8] = 0.1349112434378610;
    qr.points[8] = Barycentric3D{0.0937556561159491,
                                 0.4062443438840510,
                                 0.0937556561159491,
                                 0.4062443438840510};

    qr.weights[9] = 0.1349112434378610;
    qr.points[9] = Barycentric3D{0.0937556561159491,
                                 0.0937556561159491,
                                 0.4062443438840510,
                                 0.4062443438840510};

    return qr;
  }

  LOG_ERR("Quadrature rules of order 7 and higher have not been implemented.")
}

}
