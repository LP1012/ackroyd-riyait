#include <iostream>

#include "material.h"
#include "rectangular_region.h"
#include "simulation.h"

using namespace ar;

int main() {
  // hard-code in the regions and materials for now
  Material material_1(1, 0, 0.8, 6.4);
  Material material_2(2, 0, 0);
  Material material_3(3, 0, 0.8);

  //   // for non-scattering simulation, use:
  //   Material material_1(1, 0, 0.8, 6.4);
  //   Material material_2(2, 0, 0);
  //   Material material_3(3, 0, 0.8);

  // Code is blocked from top row to bottom row. Vacuum boundary conditions are
  // assumed for all sides.

  double target_cell_width =
      0.125 * 1.25;  // take multiples of this to get a grid convergence study

  RectangularRegion region_1(material_3, -10, 10, 5, 10, target_cell_width);

  RectangularRegion region_2(material_3, -10, -5, 1.25, 5, target_cell_width);
  RectangularRegion region_3(material_2, -5, 5, 1.25, 5, target_cell_width);
  RectangularRegion region_4(material_3, 5, 10, 1.25, 5, target_cell_width);

  RectangularRegion region_5(material_3, -10, -5, -1.25, 1.25,
                             target_cell_width);
  RectangularRegion region_6(material_2, -5, -1.25, -1.25, 1.25,
                             target_cell_width);
  RectangularRegion region_7(material_1, -1.25, 1.25, -1.25, 1.25,
                             target_cell_width);
  RectangularRegion region_8(material_2, 1.25, 5, -1.25, 1.25,
                             target_cell_width);
  RectangularRegion region_9(material_3, 5, 10, -1.25, 1.25, target_cell_width);

  RectangularRegion region_10(material_3, -10, -5, -5, -1.25,
                              target_cell_width);
  RectangularRegion region_11(material_2, -5, 5, -5, -1.25, target_cell_width);
  RectangularRegion region_12(material_3, 5, 10, -5, -1.25, target_cell_width);

  RectangularRegion region_13(material_3, -10, 10, -10, -5, target_cell_width);

  std::vector<RectangularRegion> regions = {
      region_1, region_2, region_3,  region_4,  region_5,  region_6, region_7,
      region_8, region_9, region_10, region_11, region_12, region_13};

  Simulation ar_simulation(regions, 8);
  //   ar_simulation.ExportCellsToCSV();
  ar_simulation.Run();
  ar_simulation.ExportResultsToCSV();

  return 0;
}
