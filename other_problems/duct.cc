#include <iostream>

#include "material.h"
#include "rectangular_region.h"
#include "simulation.h"

using namespace ar;

// replace main.cc with duct.cc and rename to main.cc, or copy this input
int main() {
  // hard-code in the regions and materials for now
  Material material_1(1, 0.19, 0.2, 6.4);
  Material material_2(2, 0, 0);
  Material material_3(3, 0.19, 0.2);

  //   // for non-scattering simulation, use:
  //   Material material_1(1, 0, 0.8, 6.4);
  //   Material material_2(2, 0, 0);
  //   Material material_3(3, 0, 0.8);

  // Code is blocked from top row to bottom row. Vacuum boundary conditions are
  // assumed for all sides.

  double target_cell_width =
      0.25 * 1;  // take multiples of this to get a grid convergence study

  RectangularRegion region_1(material_3, -14, -3, 0, 36, target_cell_width);
  RectangularRegion region_2(material_2, -3, 3, 0, 15, target_cell_width);
  RectangularRegion region_3(material_1, -3, 3, 15, 21, target_cell_width);
  RectangularRegion region_4(material_2, -3, 3, 21, 36, target_cell_width);
  RectangularRegion region_5(material_3, 3, 14, 0, 36, target_cell_width);

  RectangularRegion region_6(material_2, -18, -14, 0, 36, target_cell_width);
  RectangularRegion region_7(material_2, 14, 18, 0, 36, target_cell_width);

  std::vector<RectangularRegion> regions = {
      region_1, region_2, region_3, region_4, region_5, region_6, region_7};

  Simulation ar_simulation(regions, 16);
  //   ar_simulation.ExportCellsToCSV();
  ar_simulation.Run();
  ar_simulation.ExportResultsToCSV();

  return 0;
}
