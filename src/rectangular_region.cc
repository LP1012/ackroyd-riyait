#include "rectangular_region.h"

#include <cmath>
#include <ranges>

#include "point.h"

namespace ar {
RectangularRegion::RectangularRegion(const Material& material,
                                     const double xmin, const double xmax,
                                     const double ymin, const double ymax,
                                     const int nx, const int ny)
    : material_(material),
      xmin_(xmin),
      xmax_(xmax),
      ymin_(ymin),
      ymax_(ymax),
      nx_(nx),
      ny_(ny) {
  CreateCells();
}

RectangularRegion::RectangularRegion(const Material& material,
                                     const double xmin, const double xmax,
                                     const double ymin, const double ymax,
                                     const double target_delta)
    : material_(material),
      xmin_(xmin),
      xmax_(xmax),
      ymin_(ymin),
      ymax_(ymax),
      nx_(std::ceil((xmax - xmin) / target_delta)),
      ny_(std::ceil((ymax - ymin) / target_delta)) {
  CreateCells();
}

void RectangularRegion::CreateCells() {
  double dx = (xmax_ - xmin_) / static_cast<double>(nx_);
  double dy = (ymax_ - ymin_) / static_cast<double>(ny_);

  double cell_y_min = ymin_;
  for (auto j : std::views::iota(0, ny_)) {
    double cell_y_max = cell_y_min + dy;

    double cell_x_min = xmin_;
    for (auto i : std::views::iota(0, nx_)) {
      double cell_x_max = cell_x_min + dx;
      Cell new_cell = {cell_x_min, cell_x_max, cell_y_min, cell_y_max,
                       material_};
      cells_.push_back(new_cell);
      cell_x_min += dx;
    }
    cell_y_min += dy;
  }
}
}  // namespace ar
