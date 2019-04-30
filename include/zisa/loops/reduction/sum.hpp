#ifndef ZISA_SUM_HPP
#define ZISA_SUM_HPP

#include <zisa/config.hpp>
#include <zisa/loops/execution_policies.hpp>
#include <zisa/loops/range.hpp>
#include <zisa/meta/return_t.hpp>

namespace zisa::reduce {

template <class Transform>
auto sum(omp_policy, const TriangleRange &range, const Transform &transform) {
  int_t i_start = range.start_index();
  int_t i_end = range.end_index();

  using ret_type = return_t<Transform, int_t, Triangle>;
  ret_type ret = 0;

#pragma omp parallel for reduction(+ : ret)
  for (int_t i = i_start; i < i_end; ++i) {
    ret += transform(i, range.item(i));
  }

  return ret;
}

template <class Transform>
auto sum(serial_policy,
         const TriangleRange &range,
         const Transform &transform) {
  int_t i_start = range.start_index();
  int_t i_end = range.end_index();

  using ret_type = return_t<Transform, int_t, Triangle>;
  ret_type ret = 0;

  for (int_t i = i_start; i < i_end; ++i) {
    ret += transform(i, range.item(i));
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
