#ifndef ZISA_SUM_HPP
#define ZISA_SUM_HPP

#include <zisa/config.hpp>
#include <zisa/loops/execution_policies.hpp>
#include <zisa/loops/range.hpp>
#include <zisa/meta/return_t.hpp>

namespace zisa::reduce {

namespace detail {
template <class Range, class Transform>
struct return_type_trait {
  using type = decltype(std::declval<Transform>()(
      std::declval<int_t>(),
      std::declval<Range>().item(std::declval<int_t>())));
};

template <class Range, class Transform>
using return_t = typename return_type_trait<Range, Transform>::type;

}

template <class Range, class Transform>
auto sum(omp_policy, const Range &range, const Transform &transform) {
  int_t i_start = range.start_index();
  int_t i_end = range.end_index();

  using ret_type = detail::return_t<Range, Transform>;
  ret_type ret = 0;

#pragma omp parallel for reduction(+ : ret)
  for (int_t i = i_start; i < i_end; ++i) {
    ret += transform(i, range.item(i));
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
