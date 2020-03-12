#ifndef ZISA_CELL_FLAGS_HPP_UUCOO
#define ZISA_CELL_FLAGS_HPP_UUCOO

namespace zisa {

struct CellFlags {
  bool interior : 1;
  bool ghost_cell : 1;

  CellFlags() : interior(true), ghost_cell(false) {}
};

}
#endif // ZISA_GRID_FLAGS_HPP
