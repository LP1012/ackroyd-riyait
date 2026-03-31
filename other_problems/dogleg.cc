#include <iostream>

#include "material.h"
#include "rectangular_region.h"
#include "simulation.h"

using namespace ar;

// replace main.cc with dogleg.cc and rename to main.cc, or copy this input
int main() {
  // hard-code in the regions and materials for now
  Material material_1(1, 0, 0.5, 1);
  Material material_2(2, 0, 0);
  Material material_3(3, 0, 0.5);

  //   // for non-scattering simulation, use:
  //   Material material_1(1, 0, 0.8, 6.4);
  //   Material material_2(2, 0, 0);
  //   Material material_3(3, 0, 0.8);

  // Code is blocked from top row to bottom row. Vacuum boundary conditions are
  // assumed for all sides.

  double target_cell_width =
      0.25 * 1;  // take multiples of this to get a grid convergence study

  RectangularRegion region_1(material_2, -14, -9, 0, 36, target_cell_width);
  RectangularRegion region_2(material_3, -9, -6, 0, 9, target_cell_width);
  RectangularRegion region_3(material_2, -9, -3, 9, 27, target_cell_width);
  RectangularRegion region_4(material_3, -9, -6, 27, 36, target_cell_width);
  RectangularRegion region_5(material_2, -6, 6, 0, 6, target_cell_width);
  RectangularRegion region_6(material_3, -6, 6, 6, 9, target_cell_width);
  RectangularRegion region_7(material_3, -3, 3, 9, 15, target_cell_width);
  RectangularRegion region_8(material_1, -3, 3, 15, 21, target_cell_width);
  RectangularRegion region_9(material_3, -3, 3, 21, 27, target_cell_width);
  RectangularRegion region_10(material_3, -6, 6, 27, 30, target_cell_width);
  RectangularRegion region_11(material_2, -6, 6, 30, 36, target_cell_width);
  RectangularRegion region_12(material_3, 6, 9, 0, 9, target_cell_width);
  RectangularRegion region_13(material_2, 3, 9, 9, 27, target_cell_width);
  RectangularRegion region_14(material_3, 6, 9, 27, 36, target_cell_width);
  RectangularRegion region_15(material_2, 9, 14, 0, 36, target_cell_width);

  RectangularRegion region_16(material_3, -18, -14, 0, 36, target_cell_width);
  RectangularRegion region_17(material_3, 14, 18, 0, 36, target_cell_width);

  std::vector<RectangularRegion> regions = {
      region_1,  region_2,  region_3,  region_4,  region_5,  region_6,
      region_7,  region_8,  region_9,  region_10, region_11, region_12,
      region_13, region_14, region_15, region_16, region_17};

  Simulation ar_simulation(regions, 16);
  //   ar_simulation.ExportCellsToCSV();
  ar_simulation.Run();
  ar_simulation.ExportResultsToCSV();

  return 0;
}
