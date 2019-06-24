#ifndef ZISA_SYMMETRIC_CHOICES_HPP_JZEWO
#define ZISA_SYMMETRIC_CHOICES_HPP_JZEWO

#include <utility>

#include <zisa/config.hpp>

namespace zisa {

template <class Vector>
class StrictSymmetricChoices {
public:
  class EndIterator {};

  class Iterator {
  public:
    Iterator(const Vector &vector) : vector_(vector) {}

    std::pair<int_t, int_t> operator*() const {
      return {vector_[d1], vector_[d2]};
    }

    void operator++() {
      ++d2;
      if (d2 >= vector_.size()) {
        ++d1;
        d2 = d1 + 1;
      }
    }

    bool operator==(const EndIterator &) const {
      return d1 + 1 >= vector_.size();
    }
    bool operator!=(const EndIterator &other) const {
      return !((*this) == other);
    }

  private:
    int_t d1 = 0ul;
    int_t d2 = 1ul;
    const Vector &vector_;
  };

  explicit StrictSymmetricChoices(const Vector &vector) : vector_(vector) {}
  Iterator begin() const { return Iterator(vector_); }
  EndIterator end() const { return EndIterator{}; }

private:
  const Vector &vector_;
};

template <class Vector>
StrictSymmetricChoices<Vector> strict_symmetric_choices(const Vector &vector) {
  return StrictSymmetricChoices<Vector>(vector);
}

}
#endif // ZISA_SYMMETRIC_CHOICES_HPP
