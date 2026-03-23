#ifndef ACKROYD_RIYAIT_MATERIAL_H_
#define ACKROYD_RIYAIT_MATERIAL_H_

#include <cmath>

namespace ar {
class Material {
 public:
  Material(const double scattering_xs, const double total_xs,
           const double source = 0);

  double IsotropicSource() const { return source_ / 4.0 / M_PI; }
  double scattering_xs() const { return scattering_xs_; }
  double total_xs() const { return total_xs_; }

 private:
  const double scattering_xs_;
  const double total_xs_;
  const double source_;
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_MATERIAL_H_