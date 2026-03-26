#include "gauss-legendre-chebyshev.h"

#include <cmath>

#include "gauss_legendre.h"

namespace ar {

GaussLegendreChebyshev::GaussLegendreChebyshev(const int n_ordinates)
    : n_ordinates_(n_ordinates) {
  GaussLegendre gl_quad(n_ordinates_);
  std::vector<double> gl_abscissae = gl_quad.GetAbscissae();
  std::vector<double> gl_weights = gl_quad.GetWeights();

  std::vector<double> chebyshev_abscissae = ComputeChebyshevAbscissae();
  std::vector<double> chebyshev_weights =
      ComputeChebyshevWeights(chebyshev_abscissae);

  std::vector<double> etas =
      ComputeYDirectionCosines(gl_quad.GetAbscissae(), chebyshev_abscissae);

  CreateTriples(gl_abscissae, gl_weights, etas, chebyshev_weights);
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
        std::sqrt(1.0 - std::pow(gl_abscissae[i], 2)) * std::cos(omega);
    y_direction_cosines.push_back(eta);
  }
  return y_direction_cosines;
}

void GaussLegendreChebyshev::CreateTriples(
    std::vector<double>& gl_abscissae, std::vector<double>& gl_weights,
    std::vector<double>& etas, std::vector<double>& chebyshev_weights) {
  for (auto i = 0; i < n_ordinates_; i++) {
    for (auto j = 0; j < n_ordinates_; j++) {
      GLCTriplet triplet;
      double total_weight = gl_weights[i] * chebyshev_weights[j];
      triplet.weight = total_weight;
      triplet.mu = gl_abscissae[i];
      triplet.eta = etas[j];

      mu_eta_weight_.push_back(triplet);
    };
  }
}
}  // namespace ar
