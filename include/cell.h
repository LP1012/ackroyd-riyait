#ifndef ACKROYD_RIYAIT_CELL_H_
#define ACKROYD_RIYAIT_CELL_H_

#include "material.h"
#include "point.h"

namespace ar {
class Cell {
 public:
  Cell(const Point cell_center, const double dx, const double dy,
       const Material& material);

 private:
  const Point cell_center_;
  const double dx_;
  const double dy_;
  const Material& material_;
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_CELL_H_