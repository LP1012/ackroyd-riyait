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
      cell_center_(ComputeCellCenter(xmin, xmax - xmin, ymin, ymax - ymin)) {
  cell_source_ = material_.IsotropicSource();  // initialize source value
}

Point Cell::ComputeCellCenter(const double xmin, const double dx,
                              const double ymin, const double dy) {
  double x_coord = xmin + dx / 2.0;
  double y_coord = ymin + dy / 2.0;
  return Point(x_coord, y_coord);
}

void Cell::SetCellSource(const double cell_scalar_flux) {
  cell_source_ = material_.scattering_xs() * cell_scalar_flux +
                 material_.IsotropicSource();  // check this!!!
}
}  // namespace ar
