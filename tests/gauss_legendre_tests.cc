#include <gtest/gtest.h>

#include "gauss_legendre.h"

using namespace ar;

TEST(GaussLegendreTest, Constructor) { EXPECT_NO_THROW(GaussLegendre gl(2)); }