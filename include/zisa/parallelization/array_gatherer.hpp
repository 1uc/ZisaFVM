#ifndef ZISA_ARRAY_GATHERER_HPP_3D8E8
#define ZISA_ARRAY_GATHERER_HPP_3D8E8

#include <zisa/config.hpp>
#include <zisa/memory/array_view.hpp>

namespace zisa {

/// Gather a distributed array.
/** The array is, distributed along the first axis, i.e. for all k & l
 *   x(i, k)   i =  0, ...,  9 lives on rank 0
 *   x(i, k)   i = 10, ..., 17 lives on rank 1
 *  etc.
 *
 *  Task: Gather the distributed array on a subset of the total ranks.
 */
template <class T, int n_dims>
class ArrayGatherer {
protected:
  using const_view_t = array_const_view<T, n_dims, row_major>;
  using view_t = array_view<T, n_dims, row_major>;

public:
  virtual ~ArrayGatherer() = default;

  virtual void send(const const_view_t &const_view) const = 0;
  virtual void receive(const view_t &view) const = 0;
  virtual void copy_local_patch(const view_t &global,
                                const const_view_t &local) const = 0;

  /// Is (part of) the data gathered on this rank?
  virtual bool is_this_rank_gathering() const = 0;
};

}
#endif // ZISA_ARRAY_GATHERER_HPP
