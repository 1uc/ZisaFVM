// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_FOR_EACH_HPP
#define ZISA_FOR_EACH_HPP

#include <zisa/loops/execution_policies.hpp>
#include <zisa/loops/range.hpp>
#include <zisa/parallelization/omp.h>

namespace zisa {

template <class Range, class Body>
void for_each(omp_policy, const Range &range, const Body &body) {
  auto i0 = range.start_index();
  auto i_end = range.end_index();

#if ZISA_HAS_OPENMP == 1
#pragma omp parallel for ZISA_OMP_FOR_SCHEDULE_DEFAULT
#endif
  for (auto i = i0; i < i_end; ++i) {
    if constexpr (range_traits<Range>::has_item) {
      body(i, range.item(i));
    } else {
      body(i);
    }
  }
}

template <class Range, class Body>
void for_each(serial_policy, const Range &range, const Body &body) {
  auto i0 = range.start_index();
  auto i_end = range.end_index();

  for (auto i = i0; i < i_end; ++i) {
    if constexpr (range_traits<Range>::has_item) {
      body(i, range.item(i));
    } else {
      body(i);
    }
  }
}

template <class Range, class Body>
void for_each(Range &&range, Body &&body) {
  for_each(default_execution_policy{},
           std::forward<Range>(range),
           std::forward<Body>(body));
}

}

#endif // ZISA_FOR_EACH_HPP
