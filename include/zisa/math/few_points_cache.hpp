// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_FEW_POINTS_CACHE_HPP_NPQU
#define ZISA_FEW_POINTS_CACHE_HPP_NPQU

#include <zisa/config.hpp>

#include <algorithm>
#include <functional>
#include <vector>
#include <zisa/loops/reduction/min.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

template <class FX>
class FewPointsCache {
public:
  FewPointsCache() = default;
  FewPointsCache(array<XYZ, 1> points_)
      : points(std::move(points_)),
        hashes(points.size()),
        cache(points.size()) {
    auto n_points = points.size();
    array<double, 1> ps(n_points);
    bool is_good = false;

    auto random_projections = std::vector<XYZ>{
        XYZ{0.4812687223074247, 0.4901126875380488, 0.748396406012264},
        XYZ{0.5893000897592624, 0.2328024062253076, 0.9678789902593754},
        XYZ{0.312442457165149, 0.3106755425491958, 0.9489204175487449},
        XYZ{0.637790969749397, 0.0754455129073677, 0.3403685854340376},
        XYZ{0.5338831876857503, 0.9213227083629082, 0.5431566971029622},
        XYZ{0.2385209039087733, 0.913388571618219, 0.2053846497098911},
        XYZ{0.2996192337558674, 0.6389134569121155, 0.7896639672566169},
        XYZ{0.3261298736535007, 0.7045295132395185, 0.642961279569149},
        XYZ{0.0931963937441646, 0.1815748333117181, 0.7493070549814564},
        XYZ{0.5116768712958907, 0.9587701460773556, 0.6998062953160484}};

    for (const auto &rp : random_projections) {
      for (int_t i = 0; i < n_points; ++i) {
        ps[i] = zisa::dot(rp, points[i]);
      }

      std::sort(ps.begin(), ps.end());

      double dx_max = ps[n_points - 1] - ps[0];
      double dx_min
          = zisa::reduce::min(serial_policy{},
                              index_range(0, n_points - 1),
                              [&ps](int_t i) { return ps[i + 1] - ps[i]; });

      // Are the points are reasonably well separated?
      if (dx_min > dx_max / 100000.0) {
        atol = 1e-10 * dx_max;
        projection = rp;
        is_good = true;
        break;
      }
    }

    LOG_ERR_IF(!is_good, "Failed to find a suitable projection.");

    std::sort(points.begin(), points.end(), [this](const XYZ &x, const XYZ &y) {
      return zisa::dot(projection, x) < zisa::dot(projection, y);
    });

    for (int_t i = 0; i < n_points; ++i) {
      hashes[i] = zisa::dot(projection, points[i]);
    }
  }

  void update(const std::function<FX(const XYZ &x)> &f_) {
    this->f = f_;

    auto n_points = points.size();
    for (int_t i = 0; i < n_points; ++i) {
      cache[i] = f(points[i]);
    }
  }

  FX get(const XYZ &x) const {
    double hash = zisa::dot(projection, x);
    auto it = std::lower_bound(hashes.cbegin(), hashes.cend(), hash - atol);
    if (std::abs(*it - hash) <= atol) {
      return cache[integer_cast<int_t>(it - hashes.begin())];
    } else {
      return f(x);
    }
  }

private:
  std::function<FX(const XYZ &x)> f;
  array<XYZ, 1> points;
  array<double, 1> hashes;
  array<FX, 1> cache;

  XYZ projection;
  double atol;
};

}

#endif // ZISA_FEW_POINTS_CACHE_HPP
