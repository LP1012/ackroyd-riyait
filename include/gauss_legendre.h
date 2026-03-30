#ifndef ACKROYD_RIYAIT_GAUSS_LEGENDRE_H_
#define ACKROYD_RIYAIT_GAUSS_LEGENDRE_H_

#include <vector>

namespace ar {
/**
 * @brief Defines an object holding weights and abscissas for integrating on the
 * standard interval ([-1,1]) using Gauss-Legendre quadrature
 *
 */
class GaussLegendre {
 public:
  /**
   * @brief Construct a new Gauss Legendre object
   *
   * @param n_vals Number of values desired
   */
  GaussLegendre(const int n_vals);

  /**
   * @brief Get the weights
   *
   * @return std::vector<double>
   */
  const std::vector<double>& GetWeights() const { return weights_; }

  /**
   * @brief Get the abscissas
   *
   * @return std::vector<double>
   */
  const std::vector<double>& GetAbscissas() const { return abscissas_; }

 private:
  /// @brief Quadrature weights
  std::vector<double> weights_;

  /// @brief Quadrature abscissas
  std::vector<double> abscissas_;
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_GAUSS_LEGENDRE_H_