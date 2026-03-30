#include "material.h"

namespace ar {
Material::Material(const int id, const double scattering_xs,
                   const double total_xs, const double source)
    : id_(id),
      scattering_xs_(scattering_xs),
      total_xs_(total_xs),
      source_(source) {}
}  // namespace ar
