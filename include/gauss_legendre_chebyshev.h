#ifndef ACKROYD_RIYAIT_GAUSS_LEGENDRE_CHEBYSHEV_H_
#define ACKROYD_RIYAIT_GAUSS_LEGENDRE_CHEBYSHEV_H_

#include <vector>

namespace ar {

/**
 * @brief Defines a struct holding a triplet of information needed for sweeping
 * and integrating
 *
 */
struct GLCTriplet {
  double mu;      /// x-direction cosine
  double eta;     /// y-direction cosine
  double weight;  /// quadrature weight (composite)
};

/**
 * @brief Defines a quadrature object that will integrate on the unit sphere
 * using Gauss-Legendre and Gauss-Chebyshev simultaneously
 *
 */
class GaussLegendreChebyshev {
 public:
  /**
   * @brief Construct a new Gauss Legendre Chebyshev object
   *
   * @param n_ordinates Number of ordinates per quadrature set
   */
  GaussLegendreChebyshev(const int n_ordinates);

  const std::vector<GLCTriplet> &GetTriples() const { return mu_eta_weight_; }

  int n_ordinates() const { return n_ordinates_; }

 private:
  const int n_ordinates_;  /// number of ordinates per quadrature set

  std::vector<GLCTriplet>
      mu_eta_weight_;  /// vector of GLCTriplets containing quadrature info
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_GAUSS_LEGENDRE_CHEBYSHEV_H_