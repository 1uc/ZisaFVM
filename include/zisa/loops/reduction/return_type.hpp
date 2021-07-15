// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_RETURN_TYPE_HPP
#define ZISA_RETURN_TYPE_HPP

#include <zisa/config.hpp>

#include <utility>
#include <zisa/meta/return_t.hpp>

namespace zisa {
namespace reduce {

namespace detail {
template <class Range,
          class Transform,
          bool has_item = range_traits<Range>::has_item>
struct return_type_trait;

template <class Range, class Transform>
struct return_type_trait<Range, Transform, true> {
  using type = decltype(std::declval<Transform>()(
      std::declval<int_t>(),
      std::declval<Range>().item(std::declval<int_t>())));
};

template <class Range, class Transform>
struct return_type_trait<Range, Transform, false> {
  using type = decltype(std::declval<Transform>()(std::declval<int_t>()));
};

template <class Range, class Transform>
using return_t = typename return_type_trait<Range, Transform>::type;

}

}
}

#endif // ZISA_RETURN_TYPE_HPP
