#ifndef ACKROYD_RIYAIT_RECTANGULAR_REGION_H_
#define ACKROYD_RIYAIT_RECTANGULAR_REGION_H_

#include <vector>

#include "cell.h"
#include "material.h"

namespace ar {
class RectangularRegion {
 public:
  RectangularRegion(const Material& material, const double xmin,
                    const double xmax, const double ymin, const double ymax,
                    const unsigned int nx, const unsigned int ny);

 private:
  const Material& material_;
  const double xmin_;
  const double xmax_;
  const double ymin_;
  const double ymax_;
  const unsigned int nx_;
  const unsigned int ny_;

  std::vector<Cell> cells_;
  void CreateCells();
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_RECTANGULAR_REGION_H_