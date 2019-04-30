#ifndef ZISA_REDUCTION_ALL_DNWLI_HPP
#define ZISA_REDUCTION_ALL_DNWLI_HPP

#include <memory>

#include <zisa/config.hpp>
#include <zisa/loops/execution_policies.hpp>
#include <zisa/loops/range.hpp>

namespace zisa::reduce {

template <class Dereference, class Body>
bool all(omp_policy,
         const Range<PlainIndexRange, Dereference> &range,
         const Body &body) {
  auto i0 = range.start_index();
  auto i_end = range.end_index();

  bool is_good = true;

#pragma omp parallel for reduction(&& : is_good)
  for (auto i = i0; i < i_end; ++i) {
    is_good = is_good && body(i, range.item(i));
  }

  return is_good;
}

template <class Dereference, class Body>
bool all(serial_policy,
         const Range<PlainIndexRange, Dereference> &range,
         const Body &body) {

  auto i0 = range.start_index();
  auto i_end = range.end_index();

  for (auto i = i0; i < i_end; ++i) {
    if (!body(i, range.item(i))) {
      return false;
    }
  }

  return true;
}

template <class Range, class Body>
bool all(Range &&range, Body &&body) {
  return all(default_execution_policy{},
             std::forward<Range>(range),
             std::forward<Body>(body));
}

}
#endif // ZISA_ALL_HPP
