#ifndef ZISA_ARRAY_STENCIL_FAMILY_HPP_OOLQUI
#define ZISA_ARRAY_STENCIL_FAMILY_HPP_OOLQUI

#include <zisa/memory/array.hpp>
#include <zisa/reconstruction/stencil_family.hpp>

namespace zisa {

struct stencil_family_dispatch_tag {};

template <>
struct array_save_traits<StencilFamily> {
  using dispatch_tag = stencil_family_dispatch_tag;
};

void save(HierarchicalWriter &writer,
          const array_const_view<StencilFamily, 1, row_major> &stencils,
          const std::string &tag,
          stencil_family_dispatch_tag);

}
#endif // ZISA_ARRAY_STENCIL_FAMILY_HPP
