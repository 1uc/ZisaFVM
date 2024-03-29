// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_RANGE_HPP
#define ZISA_RANGE_HPP

#include <experimental/type_traits>

#include <zisa/config.hpp>

namespace zisa {

template <class IndexRange, class Dereference>
class Range : public IndexRange, public Dereference {
public:
  Range(const IndexRange &index_range, const Dereference &dereference)
      : IndexRange(index_range), Dereference(dereference) {}
};

class PlainIndexRange {
private:
  class EndIterator {};

  class Iterator {
  public:
    explicit inline Iterator(int_t i0, int_t i_end) : i(i0), i_end(i_end) {}

    inline void operator++() { i++; }
    inline int_t operator*() const { return i; }

    inline bool operator!=(const EndIterator &) const { return i < i_end; }
    inline bool operator==(const EndIterator &it) const {
      return !(*this != it);
    }

  private:
    int_t i;
    int_t i_end;
  };

public:
  explicit PlainIndexRange(int_t i_end) : i0(0), i_end(i_end) {}

  PlainIndexRange(int_t i0, int_t i_end) : i0(i0), i_end(i_end) {
    assert(i0 <= i_end);
  }

  Iterator begin() const { return Iterator(i0, i_end); }
  EndIterator end() const { return EndIterator{}; }

  int_t start_index() const { return i0; }
  int_t end_index() const { return i_end; }

private:
  int_t i0;
  int_t i_end;
};

class NoDereference {
public:
  static constexpr bool has_item() { return false; }
};

using FlatRange = Range<PlainIndexRange, NoDereference>;

template <class T>
FlatRange flat_range(const T &t) {
  return FlatRange(PlainIndexRange(0, t.size()), NoDereference{});
}

inline PlainIndexRange index_range(int_t i_start, int_t i_end) {
  return PlainIndexRange(i_start, i_end);
}

inline auto index_range(int_t n) -> decltype(auto) { return index_range(0, n); }

template <class Range>
struct range_traits {
  template <class T>
  using item_t = decltype(std::declval<T>().item(std::declval<int_t>()));

  static constexpr bool has_item
      = std::experimental::is_detected_v<item_t, Range>;
};

}

#endif // ZISA_RANGE_HPP
