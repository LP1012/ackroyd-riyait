#include "simulation.h"

namespace ar {
Simulation::Simulation(const std::vector<RectangularRegion> &regions,
                       const double si_tolerance)
    : regions_(regions), si_tolerance_(si_tolerance) {
  std::vector<Cell> flattened_cells = PullCellsFromRegions();
  cells_ = SortCells(flattened_cells);
}

std::vector<Cell> Simulation::PullCellsFromRegions() {
  std::vector<Cell> flattened_cells;
  for (auto &region : regions_) {
    for (auto &cell : region.cells()) flattened_cells.push_back(cell);
  }
  return flattened_cells;
}

std::vector<std::vector<Cell>> Simulation::SortCells(
    const std::vector<Cell> &flattened_cells) {}

double Simulation::NorthEastSweepStep(const double x_cosine,
                                      const double y_cosine, Cell &cell) {
  double numerator = cell.cell_source() +
                     2.0 * x_cosine / cell.dx() * cell.west_flux() +
                     2.0 * y_cosine / cell.dy() * cell.south_flux();
  double denominator = 2.0 * x_cosine / cell.dx() + 2.0 * y_cosine / cell.dy() +
                       cell.material().total_xs();
  return numerator / denominator;
}
}  // namespace ar
