#include "cell.h"

namespace ar {
Cell::Cell(const double xmin, const double xmax, const double ymin,
           const double ymax, const Material& material)
    : xmin_(xmin),
      xmax_(xmax),
      ymin_(ymin),
      ymax_(ymax),
      material_(material),
      dx_(xmax - xmin),
      dy_(ymax - ymin),
      cell_center_(ComputeCellCenter()) {}

Point Cell::ComputeCellCenter() {
  double x_coord = xmin_ + dx_ / 2.0;
  double y_coord = ymin_ + dy_ / 2.0;
  return Point(x_coord, y_coord);
}
}  // namespace ar
