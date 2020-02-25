#ifndef ZISA_ANY_HPP_CUPQZLj
#define ZISA_ANY_HPP_CUPQZLj

#include <memory>

#include <zisa/loops/reduction/all.hpp>

namespace zisa::reduce {

template <class Policy, class Range, class Predicate>
bool any(Policy policy, Range &&range, const Predicate &predicate) {
  return !all(policy, std::forward<Range>(range), [&predicate](const auto &x) {
    return !predicate(x);
  });
}

template <class Range, class Predicate>
bool any(Range &&range, Predicate &&predicate) {
  return any(default_execution_policy{},
             std::forward<Range>(range),
             std::forward<Predicate>(predicate));
}

}

#endif // ZISA_ANY_HPP
