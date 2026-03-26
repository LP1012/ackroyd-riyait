#ifndef ACKROYD_RIYAIT_GAUSS_LEGENDRE_CHEBYSHEV_H_
#define ACKROYD_RIYAIT_GAUSS_LEGENDRE_CHEBYSHEV_H_

#include <vector>

namespace ar {
class GaussLegendreChebyshev {
 public:
  GaussLegendreChebyshev(const int n_ordinates);

 private:
  const int n_ordinates_;

  std::vector<double> x_direction_cosines;
  std::vector<double> y_direction_cosines_;

  std::vector<double> chebyshev_abscissae_;
  std::vector<double> chebyshev_weights_;

  void SetChebyshevAbscissae(const int n_values);
  double ComputeChebyshevAbscissa(const int number);
  void SetChebyshevWeights();

  void SetEta();
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_GAUSS_LEGENDRE_CHEBYSHEV_H_