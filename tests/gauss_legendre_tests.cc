#include <gtest/gtest.h>

#include "gauss_legendre.h"

using namespace ar;

TEST(GaussLegendreTest, Constructor) { EXPECT_NO_THROW(GaussLegendre gl(2)); }

TEST(GaussLegendreTest, WeightSum) {
  for (auto i = 2; i <= 120; i += 2) {
    GaussLegendre gl(i);

    auto weights = gl.GetWeights();

    double sum = 0;
    for (auto weight : weights) sum += weight;

    EXPECT_DOUBLE_EQ(sum, 2.0);
  }
}