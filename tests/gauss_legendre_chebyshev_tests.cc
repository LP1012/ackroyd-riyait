#include <gtest/gtest.h>

#include <cmath>

#include "gauss_legendre_chebyshev.h"

using namespace ar;

TEST(GaussLegendreChebyshevTest, Constructor) {
  EXPECT_NO_THROW(GaussLegendreChebyshev gl(2));
}

TEST(GaussLegendreChebyshevTest, WeightSum) {
  for (auto i = 2; i <= 120; i += 2) {
    GaussLegendreChebyshev glc(i);

    auto& triples = glc.GetTriples();

    double sum = 0;
    for (auto& triple : triples) sum += triple.weight;

    EXPECT_NEAR(sum, 4.0 * M_PI, 1e-12);
  }
}