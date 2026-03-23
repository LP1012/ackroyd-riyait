#include "cell.h"

namespace ar {
Cell::Cell(const Point cell_center, const double dx, const double dy,
           const Material& material)
    : cell_center_(cell_center), dx_(dx), dy_(dy), material_(material) {}
}  // namespace ar
