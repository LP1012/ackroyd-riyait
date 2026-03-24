#include "simulation.h"

#include <cmath>

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

double Simulation::SweepStep(const double x_cosine, const double y_cosine,
                             Cell &cell) {
  double east_west_flux = cell.west_flux();
  double north_south_flux = cell.south_flux();
  if (x_cosine < 0) east_west_flux = cell.east_flux();

  if (y_cosine < 0) north_south_flux = cell.north_flux();

  double numerator = cell.cell_source() +
                     2.0 * std::abs(x_cosine) / cell.dx() * east_west_flux +
                     2.0 * std::abs(y_cosine) / cell.dy() * north_south_flux;
  double denominator = 2.0 * std::abs(x_cosine) / cell.dx() +
                       2.0 * std::abs(y_cosine) / cell.dy() +
                       cell.material().total_xs();
  return numerator / denominator;
}

}  // namespace ar
