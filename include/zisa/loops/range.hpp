#ifndef ZISA_RANGE_HPP
#define ZISA_RANGE_HPP

#include <zisa/config.hpp>

namespace zisa {

template <class IndexRange, class Dereference>
class Range : public IndexRange, public Dereference {
public:
  Range(const IndexRange &index_range, const Dereference &dereference)
      : IndexRange(index_range), Dereference(dereference) {}
};

class PlainIndexRange {
public:
  PlainIndexRange(int_t i0, int_t i_end) : i0(i0), i_end(i_end) {
    assert(i0 <= i_end);
  }

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

}

#endif // ZISA_RANGE_HPP
