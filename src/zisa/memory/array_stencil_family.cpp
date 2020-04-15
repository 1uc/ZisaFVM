#include <zisa/memory/array_stencil_family.hpp>

namespace zisa {

void save(HDF5Writer &writer,
          const array<StencilFamily, 1, row_major> &stencils,
          const std::string &tag,
          stencil_family_dispatch_tag) {

  bool needs_final_close = false;
  if (tag != "") {
    writer.open_group(tag);
    needs_final_close = true;
  }

  auto n_cells = stencils.shape(0);
  for (int_t i = 0; i < n_cells; ++i) {
    writer.open_group(string_format("%d", i));

    const auto &sf = stencils[i];
    for (int_t k = 0; k < sf.size(); ++k) {
      save(writer, sf[k].global(), string_format("%d", k));
    }
    writer.close_group();
  }

  if (needs_final_close) {
    writer.close_group();
  }
}

}