#include "rectangular_region.h"

namespace ar {
RectangularRegion::RectangularRegion(const Material& material,
                                     const double xmin, const double xmax,
                                     const double ymin, const double ymax,
                                     const unsigned int nx,
                                     const unsigned int ny)
    : material_(material),
      xmin_(xmin),
      xmax_(xmax),
      ymin_(ymin),
      ymax_(ymax),
      nx_(nx),
      ny_(ny) {
  CreateCells();
}
}  // namespace ar
