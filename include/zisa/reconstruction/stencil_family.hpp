#ifndef STENCIL_FAMILY_H_4HHJ0
#define STENCIL_FAMILY_H_4HHJ0

#include <zisa/config.hpp>
#include <zisa/reconstruction/stencil.hpp>
#include <zisa/reconstruction/stencil_family_params.hpp>

namespace zisa {

class StencilFamily {
private:
  array<Stencil, 1> stencils_;

public:
  StencilFamily() = default;
  StencilFamily(const Grid &grid,
                int_t i_cell,
                const StencilFamilyParams &params);

  /// Returns the k-th stencil.
  const Stencil &operator[](int_t k) const;

  /// Mapping from local to global indices.
  const std::vector<int_t> &local2global() const;

  /// Return the number of stencils.
  int_t size() const;

  /// Return the size of the combined stencil.
  int_t combined_stencil_size() const;

  /// Return the index of the highest order central stencil.
  int_t highest_order_stencil() const;

  /// Returns the highest order possible for the given stencils.
  int order() const;

  auto begin() -> decltype(stencils_.begin());
  auto begin() const -> decltype(stencils_.begin());

  auto end() -> decltype(stencils_.end());
  auto end() const -> decltype(stencils_.end());

  template<class F>
  void apply_permutation(const F &f) {
    for(auto &s : (*this)) {
      s.apply_permutation(f);
    }
    for(auto &i : l2g) {
      i = f(i);
    }
  }

protected:
  /// This sets the order of every stencil to 1.
  void truncate_all_stencils_to_first_order(int_t i_cell);

private:
  std::vector<int_t> l2g;

  int_t combined_stencil_size_;
  int order_;
  int_t k_high_;
};

bool operator==(const StencilFamily &lhs, const StencilFamily &rhs);
bool operator!=(const StencilFamily &lhs, const StencilFamily &rhs);

array<StencilFamily, 1>
compute_stencil_families(const Grid &grid,
                         const StencilFamilyParams &params);

} // namespace zisa
#endif /* end of include guard */
