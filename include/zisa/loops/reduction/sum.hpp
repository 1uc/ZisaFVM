// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_SUM_HPP
#define ZISA_SUM_HPP

#include <zisa/config.hpp>
#include <zisa/loops/execution_policies.hpp>
#include <zisa/loops/range.hpp>
#include <zisa/loops/reduction/return_type.hpp>

namespace zisa::reduce {

template <class Range, class Transform>
auto sum(omp_policy, const Range &range, const Transform &transform) {
  int_t i_start = range.start_index();
  int_t i_end = range.end_index();

  using ret_type = detail::return_t<Range, Transform>;
  ret_type ret = 0;

#if ZISA_HAS_OPENMP == 1
#pragma omp parallel for reduction(+ : ret)
#endif
  for (int_t i = i_start; i < i_end; ++i) {
    if constexpr (range_traits<Range>::has_item) {
      ret += transform(i, range.item(i));
    } else {
      ret += transform(i);
    }
  }

  return ret;
}

template <class Range, class Transform>
auto sum(serial_policy, const Range &range, const Transform &transform) {
  int_t i_start = range.start_index();
  int_t i_end = range.end_index();

  using ret_type = detail::return_t<Range, Transform>;
  ret_type ret = 0;

  for (int_t i = i_start; i < i_end; ++i) {
    if constexpr (range_traits<Range>::has_item) {
      ret += transform(i, range.item(i));
    } else {
      ret += transform(i);
    }
  }

  return ret;
}

template <class Range, class Transform>
auto sum(Range &&range, Transform &&transform) -> decltype(auto) {

  return zisa::reduce::sum(default_execution_policy{},
                           std::forward<Range>(range),
                           std::forward<Transform>(transform));
}

}

#endif // ZISA_SUM_HPP
