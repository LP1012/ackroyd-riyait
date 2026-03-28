#ifndef ACKROYD_RIYAIT_CELL_H_
#define ACKROYD_RIYAIT_CELL_H_

#include "material.h"
#include "point.h"

namespace ar {
class Cell {
 public:
  Cell(const double xmin, const double xmax, const double ymin,
       const double ymax, const Material& material);

  double cell_source() const { return cell_source_; }
  double west_flux() const { return west_flux_; }
  double north_flux() const { return north_flux_; }
  double south_flux() const { return south_flux_; }
  double east_flux() const { return east_flux_; }
  double dx() const { return dx_; }
  double dy() const { return dy_; }
  Point cell_center() const { return cell_center_; }
  Material material() const { return material_; }

  void SetEastFlux() { east_flux_ = 2 * center_flux_ - west_flux_; }
  void SetWestFlux() { west_flux_ = 2 * center_flux_ - east_flux_; }
  void SetNorthFlux() { north_flux_ = 2 * center_flux_ - south_flux_; }
  void SetSouthFlux() { south_flux_ = 2 * center_flux_ - north_flux_; }
  void SetCenterFlux(const double center_flux) { center_flux_ = center_flux; }

 private:
  const Point cell_center_;
  const double xmin_;
  const double xmax_;
  const double ymin_;
  const double ymax_;
  const double dx_;
  const double dy_;
  const Material& material_;

  Point ComputeCellCenter(const double xmin, const double dx, const double ymin,
                          const double dy);

  double west_flux_ = 0;
  double north_flux_ = 0;
  double south_flux_ = 0;
  double east_flux_ = 0;
  double center_flux_ = 0;
  double cell_source_ = 0;

  void SetCellSource(const double cell_scalar_flux);
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_CELL_H_
