#ifndef ZISA_REDUCTION_ALL_DNWLI_HPP
#define ZISA_REDUCTION_ALL_DNWLI_HPP

#include <memory>

#include <zisa/config.hpp>
#include <zisa/loops/execution_policies.hpp>
#include <zisa/loops/range.hpp>

namespace zisa::reduce {

template <class Range, class Predicate>
bool all(omp_policy, const Range &range, const Predicate &predicate) {
  auto i0 = range.start_index();
  auto i_end = range.end_index();

  bool is_good = true;

#pragma omp parallel for reduction(&& : is_good)
  for (auto i = i0; i < i_end; ++i) {
    is_good = is_good && predicate(i, range.item(i));

    if constexpr (range_traits<Range>::has_item) {
      is_good = is_good && predicate(i, range.item(i));
    } else {
      is_good = is_good && predicate(i);
    }
  }

  return is_good;
}

template <class Range, class Predicate>
bool all(serial_policy, const Range &range, const Predicate &predicate) {

  auto i0 = range.start_index();
  auto i_end = range.end_index();

  for (auto i = i0; i < i_end; ++i) {
    if constexpr (range_traits<Range>::has_item) {
      if (!predicate(i, range.item(i))) {
        return false;
      }
    } else {
      if (!predicate(i)) {
        return false;
      }
    }
  }

  return true;
}

template <class Range, class Predicate>
bool all(Range &&range, Predicate &&predicate) {
  return all(default_execution_policy{},
             std::forward<Range>(range),
             std::forward<Predicate>(predicate));
}

}
#endif // ZISA_ALL_HPP
