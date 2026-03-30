#include "gauss_legendre_chebyshev.h"

#include <cmath>

#include "gauss_legendre.h"

namespace ar {

GaussLegendreChebyshev::GaussLegendreChebyshev(const int n_ordinates)
    : n_ordinates_(n_ordinates) {
  GaussLegendre gl_quad(n_ordinates_);
  std::vector<double> gl_abscissae = gl_quad.GetAbscissas();
  std::vector<double> gl_weights = gl_quad.GetWeights();

  const double chebyshev_weight =
      2.0 * M_PI / static_cast<double>(n_ordinates_);

  for (auto i = 0; i < n_ordinates_; i++) {
    for (auto j = 0; j < n_ordinates_; j++) {
      double mu = gl_abscissae[i];
      double phi = (2.0 * j + 1.0) * M_PI / static_cast<double>(n_ordinates_);
      double eta = std::sqrt(1.0 - mu * mu) * std::sin(phi);

      GLCTriplet triplet;
      triplet.mu = mu;
      triplet.eta = eta;
      triplet.weight = gl_weights[i] * chebyshev_weight;
      mu_eta_weight_.push_back(triplet);
    }
  }
}
}  // namespace ar
