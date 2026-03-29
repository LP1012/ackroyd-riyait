#include <gtest/gtest.h>

#include "simulation.h"

using namespace ar;

TEST(SimulationTest, Constructor) {
  Material material_1(1, 1, 2);
  Material material_2(2, 0, 1, 1);
  RectangularRegion region_1(material_1, 0, 1, 0, 1, 2, 2);
  RectangularRegion region_2(material_2, 1, 2, 0, 1, 2, 2);

  std::vector<RectangularRegion> regions = {region_1, region_2};
  EXPECT_NO_THROW(Simulation simulation(regions));
}
