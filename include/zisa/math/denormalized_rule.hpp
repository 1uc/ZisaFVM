// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_DENORMALIZEDRULE_HPP_ICPDOI
#define ZISA_DENORMALIZEDRULE_HPP_ICPDOI

#include <zisa/config.hpp>
#include <zisa/io/format_as_list.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

struct DenormalizedRule {
  array<double, 1> weights;
  array<XYZ, 1> points;
  double volume;

  DenormalizedRule() = default;
  explicit DenormalizedRule(int_t n_points);
  DenormalizedRule(const DenormalizedRule &qr) = default;
  DenormalizedRule(DenormalizedRule &&qr) = default;

  DenormalizedRule &operator=(const DenormalizedRule &qr) = default;
  DenormalizedRule &operator=(DenormalizedRule &&qr) = default;
};

std::string str(const DenormalizedRule &a);

bool operator==(const DenormalizedRule &a, const DenormalizedRule &b);

bool operator!=(const DenormalizedRule &a, const DenormalizedRule &b);

inline double volume(const DenormalizedRule &qr) { return qr.volume; }

template <class QR, class Domain>
DenormalizedRule denormalize(const QR &qr_hat, const Domain &domain) {
  assert(qr_hat.weights.size() == qr_hat.points.size());

  auto n_points = qr_hat.weights.size();
  auto qr = DenormalizedRule(n_points);

  auto vol = volume(domain);
  for (int_t i = 0; i < n_points; ++i) {
    qr.weights[i] = vol * qr_hat.weights[i];
    qr.points[i] = coord(domain, qr_hat.points[i]);
  }

  qr.volume = vol;
  return qr;
}

}

#endif // ZISA_DENORMALIZEDRULE_HPP
