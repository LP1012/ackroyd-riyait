#include "gauss-legendre-chebyshev.h"

#include <cmath>

namespace ar {

GaussLegendreChebyshev::GaussLegendreChebyshev(const int n_ordinates)
    : n_ordinates_(n_ordinates) {
  SetChebyshevAbscissae();
  SetChebyshevWeights();
  SetEta();
}

void GaussLegendreChebyshev::SetChebyshevAbscissae() {
  for (auto i = 0; i < n_ordinates_; i++)
    chebyshev_abscissae_.push_back(ComputeChebyshevAbscissa(i));
}

double GaussLegendreChebyshev::ComputeChebyshevAbscissa(const int number) {
  return std::cos((2.0 * static_cast<double>(number) - 1.0) /
                  (2.0 * static_cast<double>(n_ordinates_)) * M_PI);
}

void GaussLegendreChebyshev::SetChebyshevWeights() {
  double front_coeff = M_PI * M_PI / static_cast<double>(n_ordinates_);
  for (auto& value : chebyshev_abscissae_)
    chebyshev_weights_.push_back(front_coeff * std::sqrt(1.0 - value * value));
}

void GaussLegendreChebyshev::SetEta() {
  for (auto i = 0; i < n_ordinates_; i++) {
    double omega = M_PI * (chebyshev_abscissae_[i] + 1.0);
    double eta =
        std::sqrt(1.0 - std::pow(x_direction_cosines_[i], 2)) * std::cos(omega);
    y_direction_cosines_.push_back(eta);
  }
}

}  // namespace ar
