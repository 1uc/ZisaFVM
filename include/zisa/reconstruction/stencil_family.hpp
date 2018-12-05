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
  StencilFamily(const std::shared_ptr<Grid> &grid,
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

  auto begin() const -> decltype(stencils_.begin());
  auto end() const -> decltype(stencils_.end());

private:
  std::vector<int_t> l2g;

  int_t combined_stencil_size_;
  int order_;
  int_t k_high_;
};

bool operator==(const StencilFamily &lhs, const StencilFamily &rhs);
bool operator!=(const StencilFamily &lhs, const StencilFamily &rhs);

} // namespace zisa
#endif /* end of include guard */
