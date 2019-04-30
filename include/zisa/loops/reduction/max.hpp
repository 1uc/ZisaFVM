#ifndef ZISA_MAX_HPP
#define ZISA_MAX_HPP

#include <memory>

#include <zisa/config.hpp>
#include <zisa/loops/execution_policies.hpp>
#include <zisa/loops/range.hpp>
#include <zisa/meta/return_t.hpp>

namespace zisa::reduce {

template <class Transform>
auto max(omp_policy, const TriangleRange &range, const Transform &transform) {
  int_t i_start = range.start_index();
  int_t i_end = range.end_index();

  using ret_type = return_t<Transform, int_t, Triangle>;
  ret_type ret = std::numeric_limits<ret_type>::lowest();

#pragma omp parallel for reduction(max : ret)
  for (int_t i = i_start; i < i_end; ++i) {
    ret = zisa::max(ret, transform(i, range.item(i)));
  }

  return ret;
}

template <class Transform>
auto max(serial_policy,
         const TriangleRange &range,
         const Transform &transform) {
  int_t i_start = range.start_index();
  int_t i_end = range.end_index();

  using ret_type = return_t<Transform, int_t, Triangle>;
  ret_type ret = std::numeric_limits<ret_type>::lowest();

  for (int_t i = i_start; i < i_end; ++i) {
    ret = zisa::max(ret, transform(i, range.item(i)));
  }

  return ret;
}

template <class Range, class Transform>
auto max(Range &&range, Transform &&transform) -> decltype(auto) {

  return zisa::reduce::max(default_execution_policy{},
                           std::forward<Range>(range),
                           std::forward<Transform>(transform));
}

}
#endif // ZISA_MAX_HPP
