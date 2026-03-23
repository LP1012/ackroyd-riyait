#include "cell.h"

namespace ar {
Cell::Cell(const Point cell_center, const double dx, const double dy)
    : cell_center_(cell_center), dx_(dx), dy_(dy) {}
}  // namespace ar
