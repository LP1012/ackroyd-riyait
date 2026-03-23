#include "material.h"

namespace ar {
Material::Material(const double scattering_xs, const double total_xs,
                   const double source)
    : scattering_xs_(scattering_xs_), total_xs_(total_xs), source_(source) {}
}  // namespace ar
