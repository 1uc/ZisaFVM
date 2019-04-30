#ifndef ZISA_MIN_KDOW_HPP
#define ZISA_MIN_KDOW_HPP

#include <zisa/config.hpp>
#include <zisa/loops/execution_policies.hpp>
#include <zisa/loops/range.hpp>
#include <zisa/meta/return_t.hpp>

namespace zisa::reduce {

template <class Transform>
auto min(omp_policy, const TriangleRange &range, const Transform &transform) {
  int_t i_start = range.start_index();
  int_t i_end = range.end_index();

  using ret_type = return_t<Transform, int_t, Triangle>;
  ret_type ret = std::numeric_limits<ret_type>::max();

#pragma omp parallel for reduction(min : ret)
  for (int_t i = i_start; i < i_end; ++i) {
    ret = zisa::min(ret, transform(i, range.item(i)));
  }

  return ret;
}

template <class Transform>
auto min(serial_policy,
         const TriangleRange &range,
         const Transform &transform) {
  int_t i_start = range.start_index();
  int_t i_end = range.end_index();

  using ret_type = return_t<Transform, int_t, Triangle>;
  ret_type ret = std::numeric_limits<ret_type>::max();

  for (int_t i = i_start; i < i_end; ++i) {
    ret = zisa::min(ret, transform(i, range.item(i)));
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
