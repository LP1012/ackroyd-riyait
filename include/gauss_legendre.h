#ifndef ACKROYD_RIYAIT_GAUSS_LEGENDRE_H_
#define ACKROYD_RIYAIT_GAUSS_LEGENDRE_H_

#include <vector>

namespace ar {
class GaussLegendre {
 public:
  GaussLegendre(const int n_vals);

  std::vector<double> GetWeights() const { return weights_; }
  std::vector<double> GetAbscissae() const { return abscissae_; }

 private:
  std::vector<double> weights_;
  std::vector<double> abscissae_;
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_GAUSS_LEGENDRE_H_