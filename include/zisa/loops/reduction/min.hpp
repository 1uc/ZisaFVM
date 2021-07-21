// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_MIN_KDOW_HPP
#define ZISA_MIN_KDOW_HPP

#include <zisa/config.hpp>
#include <zisa/loops/execution_policies.hpp>
#include <zisa/loops/range.hpp>
#include <zisa/loops/reduction/return_type.hpp>
#include <zisa/math/comparison.hpp>

namespace zisa::reduce {

template <class Range, class Transform>
auto min(omp_policy, const Range &range, const Transform &transform) {
  int_t i_start = range.start_index();
  int_t i_end = range.end_index();

  using ret_type = detail::return_t<Range, Transform>;
  ret_type ret = std::numeric_limits<ret_type>::max();

#if ZISA_HAS_OPENMP == 1
#pragma omp parallel for reduction(min : ret)
#endif
  for (int_t i = i_start; i < i_end; ++i) {
    if constexpr (range_traits<Range>::has_item) {
      ret = zisa::min(ret, transform(i, range.item(i)));
    } else {
      ret = zisa::min(ret, transform(i));
    }
  }

  return ret;
}

template <class Range, class Transform>
auto min(serial_policy, const Range &range, const Transform &transform) {
  int_t i_start = range.start_index();
  int_t i_end = range.end_index();

  using ret_type = detail::return_t<Range, Transform>;
  ret_type ret = std::numeric_limits<ret_type>::max();

  for (int_t i = i_start; i < i_end; ++i) {
    if constexpr (range_traits<Range>::has_item) {
      ret = zisa::min(ret, transform(i, range.item(i)));
    } else {
      ret = zisa::min(ret, transform(i));
    }
  }

  return ret;
}

template <class Range, class Transform>
auto min(Range &&range, Transform &&transform) -> decltype(auto) {

  return zisa::reduce::min(default_execution_policy{},
                           std::forward<Range>(range),
                           std::forward<Transform>(transform));
}

}

#endif // ZISA_MIN_HPP
