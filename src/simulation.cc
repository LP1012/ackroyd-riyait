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
  double x_sign = x_cosine > 0 ? 1.0 : -1.0;
  double y_sign = y_cosine > 0 ? 1.0 : -1.0;
  double numerator = cell.cell_source() +
                     x_sign * 2.0 * x_cosine / cell.dx() * cell.west_flux() +
                     y_sign * 2.0 * y_cosine / cell.dy() * cell.south_flux();
  double denominator = x_sign * 2.0 * x_cosine / cell.dx() +
                       y_sign * 2.0 * y_cosine / cell.dy() +
                       cell.material().total_xs();
  return numerator / denominator;
}

double Simulation::NorthWestSweepStep(const double x_cosine,
                                      const double y_cosine, Cell &cell) {
  double x_sign = x_cosine > 0 ? 1.0 : -1.0;
  double y_sign = y_cosine > 0 ? 1.0 : -1.0;
  double numerator = cell.cell_source() +
                     x_sign * 2.0 * x_cosine / cell.dx() * cell.east_flux() +
                     y_sign * 2.0 * y_cosine / cell.dy() * cell.south_flux();
  double denominator = x_sign * 2.0 * x_cosine / cell.dx() +
                       y_sign * 2.0 * y_cosine / cell.dy() +
                       cell.material().total_xs();
  return numerator / denominator;
}
}  // namespace ar
