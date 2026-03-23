#ifndef ACKROYD_RIYAIT_CELL_H_
#define ACKROYD_RIYAIT_CELL_H_

#include "point.h"

namespace ar {
class Cell {
 public:
  Cell(const Point cell_center, const double dx, const double dy);

 private:
  const Point cell_center_;
  const double dx_;
  const double dy_;
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_CELL_H_