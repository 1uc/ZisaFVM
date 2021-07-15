// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_PERMUTATION_HPP_CAIWH
#define ZISA_PERMUTATION_HPP_CAIWH

#include <vector>

#include <zisa/config.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/memory/array_view.hpp>

namespace zisa {

struct Cycle {
  Cycle() = default;
  Cycle(const Cycle &) = default;
  Cycle(Cycle &&) = default;

  Cycle &operator=(const Cycle &) = default;
  Cycle &operator=(Cycle &&) = default;

  explicit Cycle(int_t cycle_length) : cycle(shape_t<1>(cycle_length)) {}

  int_t operator()(int_t i) const { return cycle(i); }
  int_t &operator()(int_t i) { return cycle(i); }

  int_t operator[](int_t i) const { return cycle(i); }
  int_t &operator[](int_t i) { return cycle(i); }

  int_t size() const { return cycle.size(); }

private:
  array<int_t, 1> cycle;
};

struct Permutation {
  array<Cycle, 1> cycles;
};

/// Factor the permutation into disjoint cycles.
Permutation factor_permutation(const array_const_view<int_t, 1> &sigma);

/// Apply the permutation to the array.
/** This inplace transform applies the permutation to the first axis, such that
 *
 *    1D: x_new(i) == x_old(sigma(i))
 *    2D: x_new(i, k) == x_old(sigma(i), k) for all k.
 *    3D: ... (needs to be implemented).
 *
 *  Individual elements are moved.
 */
template <class T, int n_dims>
void apply_permutation(array_view<T, n_dims> data, const Permutation &sigma);

/// Reverse the permutation of the array.
/** This inplace transform applies the reverse permutation to the first axis,
 *  such that
 *
 *    1D: x_new(sigma(i)) == x_old(i)
 *    2D: x_new(sigma(i), k) == x_old(i, k) for all k.
 *    3D: ... (needs to be implemented).
 *
 *  Individual elements are moved.
 */
template <class T, int n_dims>
void reverse_permutation(array_view<T, n_dims> data, const Permutation &sigma);



// -- Implementation -----------------------------------------------------------
template <class T>
void apply_cycle(array_view<T, 1> &data, const Cycle &cycle) {
  T first = std::move(data[cycle[0]]);

  for (int_t i = 0; i < cycle.size() - 1; ++i) {
    data[cycle[i]] = std::move(data[cycle[i + 1]]);
  }
  data[cycle[cycle.size() - 1]] = std::move(first);
}

template <class T>
void apply_cycle(array_view<T, 2> &data, const Cycle &cycle) {
  auto cycle_length = cycle.size();
  auto n_vars = data.shape(1);

  array<T, 1> first(n_vars);
  for (int_t k = 0; k < n_vars; ++k) {
    first(k) = std::move(data(cycle(0), k));
  }

  for (int_t i = 0; i < cycle_length - 1; ++i) {
    for (int_t k = 0; k < n_vars; ++k) {
      data(cycle(i), k) = std::move(data(cycle(i + 1), k));
    }
  }

  for (int_t k = 0; k < n_vars; ++k) {
    data(cycle(cycle_length - 1), k) = std::move(first(k));
  }
}

template <class T, int n_dims>
void apply_permutation(array_view<T, n_dims> data, const Permutation &sigma) {
  if(data.shape(0) == 0) {
    return;
  }

  for (const auto &cycle : sigma.cycles) {
    apply_cycle(data, cycle);
  }
}

// -- Reverse ------------------------------------------------------------------

template <class T>
void reverse_cycle(array_view<T, 1> &data, const Cycle &cycle) {
  auto cycle_size = cycle.size();
  T last = std::move(data[cycle[cycle_size-1]]);

  for (int_t i = cycle_size-1; i > 0; --i) {
    data[cycle[i]] = std::move(data[cycle[i-1]]);
  }
  data[cycle[0]] = std::move(last);
}

template <class T>
void reverse_cycle(array_view<T, 2> &data, const Cycle &cycle) {
  auto cycle_length = cycle.size();
  auto n_vars = data.shape(1);

  array<T, 1> last(n_vars);
  for (int_t k = 0; k < n_vars; ++k) {
    last(k) = std::move(data(cycle(cycle_length-1), k));
  }

  for (int_t i = cycle_length-1; i > 0; --i) {
    for (int_t k = 0; k < n_vars; ++k) {
      data(cycle(i), k) = std::move(data(cycle(i-1), k));
    }
  }

  for (int_t k = 0; k < n_vars; ++k) {
    data(cycle(0), k) = std::move(last(k));
  }
}

template <class T, int n_dims>
void reverse_permutation(array_view<T, n_dims> data, const Permutation &sigma) {
  if(data.shape(0) == 0) {
    return;
  }

  for (const auto &cycle : sigma.cycles) {
    reverse_cycle(data, cycle);
  }
}

}

#endif // ZISA_PERMUTATION_HPP
