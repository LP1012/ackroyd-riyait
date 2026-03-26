#ifndef ACKROYD_RIYAIT_GAUSS_LEGENDRE_CHEBYSHEV_H_
#define ACKROYD_RIYAIT_GAUSS_LEGENDRE_CHEBYSHEV_H_

#include <tuple>
#include <vector>

namespace ar {
class GaussLegendreChebyshev {
 public:
  GaussLegendreChebyshev(const int n_ordinates);

 private:
  const int n_ordinates_;

  std::vector<std::tuple<double, double, double>> mu_eta_weight_;

  std::vector<double> x_direction_cosines_;

  std::vector<double> ComputeChebyshevAbscissae();
  double ComputeChebyshevAbscissaValue(const int number);
  std::vector<double> ComputeChebyshevWeights(
      const std::vector<double> &chebyshev_abscissae);

  std::vector<double> ComputeYDirectionCosines(
      const std::vector<double> &gl_abscissae,
      const std::vector<double> &chebyshev_abscissae);
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_GAUSS_LEGENDRE_CHEBYSHEV_H_