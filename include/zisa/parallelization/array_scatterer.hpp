#ifndef ZISA_ARRAY_SCATTERER_HPP_ENXZIU
#define ZISA_ARRAY_SCATTERER_HPP_ENXZIU

#include <zisa/config.hpp>

#include <zisa/memory/array_view.hpp>

namespace zisa {

/// Scatters a distributed array.
/** This reverses the effects of `ArrayScatterer`. */
template <class T, int n_dims>
class ArrayScatterer {
protected:
  using const_view_t = array_const_view<T, n_dims, row_major>;
  using view_t = array_view<T, n_dims, row_major>;

public:
  virtual ~ArrayScatterer() = default;

  virtual void send(const const_view_t &const_view) const = 0;
  virtual void receive(const view_t &view) const = 0;
  virtual void copy_local_patch(const view_t &view,
                                const const_view_t &local) const = 0;

  /// Is (part of) the data scattered from this rank?
  virtual bool is_this_rank_scattering() const = 0;
};

}

#endif // ZISA_ARRAY_SCATTERER_HPP
