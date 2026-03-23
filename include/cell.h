#ifndef ACKROYD_RIYAIT_CELL_H_
#define ACKROYD_RIYAIT_CELL_H_

#include "material.h"
#include "point.h"

namespace ar {
class Cell {
 public:
  Cell(const double xmin, const double xmax, const double ymin,
       const double ymax, const Material& material);

 private:
  const Point cell_center_;
  const double xmin_;
  const double xmax_;
  const double ymin_;
  const double ymax_;
  const double dx_;
  const double dy_;
  const Material& material_;

  Point ComputeCellCenter();
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_CELL_H_