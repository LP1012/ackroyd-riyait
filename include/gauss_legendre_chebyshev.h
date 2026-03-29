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

 private:
  const int n_ordinates_;  /// number of ordinates per quadrature set

  std::vector<GLCTriplet>
      mu_eta_weight_;  /// vector of GLCTriplets containing quadrature info

  /**
   * @brief Computes the Gauss-Chebyshev abscissae on 0 to 2pi
   *
   * @return std::vector<double>
   */
  std::vector<double> ComputeChebyshevAbscissae();

  /**
   * @brief Computes single Gauss-Chebyshev abscissa value given a number of the
   * total to be generated
   *
   * @param number
   * @return double
   */
  double ComputeChebyshevAbscissaValue(const int number);

  /**
   * @brief Computes Gauss-Chebyshev weights
   *
   * @param chebyshev_abscissae
   * @return std::vector<double>
   */
  std::vector<double> ComputeChebyshevWeights(
      const std::vector<double> &chebyshev_abscissae);

  /**
   * @brief Computes y-direction cosine values ("eta" values)
   *
   * @param gl_abscissae Gauss-Legendre abscissae
   * @param chebyshev_abscissae Gauss-Chebyshev abscissae
   * @return std::vector<double>
   */
  std::vector<double> ComputeYDirectionCosines(
      const std::vector<double> &gl_abscissae,
      const std::vector<double> &chebyshev_abscissae);

  /**
   * @brief Create a GLCTriplets
   *
   * @param gl_abscissae Gauss-Legendre abscissae
   * @param gl_weights Gauss-Legendre weights
   * @param etas y-direction cosine values
   * @param chebyshev_weights Gauss-Chebyshev weights
   */
  void CreateTriples(std::vector<double> &gl_abscissae,
                     std::vector<double> &gl_weights, std::vector<double> &etas,
                     std::vector<double> &chebyshev_weights);
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_GAUSS_LEGENDRE_CHEBYSHEV_H_