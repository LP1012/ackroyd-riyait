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
                    const int nx, const int ny);
  RectangularRegion(const Material& material, const double xmin,
                    const double xmax, const double ymin, const double ymax,
                    const double target_delta);

  std::vector<Cell> cells() const { return cells_; }

 private:
  const Material& material_;
  const double xmin_;
  const double xmax_;
  const double ymin_;
  const double ymax_;
  const int nx_;
  const int ny_;

  std::vector<Cell> cells_;
  void CreateCells();
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_RECTANGULAR_REGION_H_