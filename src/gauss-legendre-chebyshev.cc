#include "gauss-legendre-chebyshev.h"

#include <cmath>

namespace ar {

GaussLegendreChebyshev::GaussLegendreChebyshev(const int n_ordinates)
    : n_ordinates_(n_ordinates) {
  std::vector<double> chebyshev_abscissae = ComputeChebyshevAbscissae();
  std::vector<double> chebyshev_weights =
      ComputeChebyshevWeights(chebyshev_abscissae);
  std::vector<double> etas =
      ComputeYDirectionCosines(gl_abscissae, chebyshev_abscissae);
}

std::vector<double> GaussLegendreChebyshev::ComputeChebyshevAbscissae() {
  std::vector<double> chebyshev_abscissae;
  for (auto i = 0; i < n_ordinates_; i++)
    chebyshev_abscissae.push_back(ComputeChebyshevAbscissaValue(i));
  return chebyshev_abscissae;
}

double GaussLegendreChebyshev::ComputeChebyshevAbscissaValue(const int number) {
  return std::cos((2.0 * static_cast<double>(number) - 1.0) /
                  (2.0 * static_cast<double>(n_ordinates_)) * M_PI);
}

std::vector<double> GaussLegendreChebyshev::ComputeChebyshevWeights(
    const std::vector<double>& chebyshev_abscissae) {
  std::vector<double> chebyshev_weights;
  double front_coeff = M_PI * M_PI / static_cast<double>(n_ordinates_);
  for (auto& value : chebyshev_abscissae)
    chebyshev_weights.push_back(front_coeff * std::sqrt(1.0 - value * value));
  return chebyshev_weights;
}

std::vector<double> GaussLegendreChebyshev::ComputeYDirectionCosines(
    const std::vector<double>& gl_abscissae,
    const std::vector<double>& chebyshev_abscissae) {
  std::vector<double> y_direction_cosines;
  for (auto i = 0; i < n_ordinates_; i++) {
    double omega = M_PI * (chebyshev_abscissae[i] + 1.0);
    double eta =
        std::sqrt(1.0 - std::pow(x_direction_cosines_[i], 2)) * std::cos(omega);
    y_direction_cosines.push_back(eta);
  }
  return y_direction_cosines;
}

}  // namespace ar
