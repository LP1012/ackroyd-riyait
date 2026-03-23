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
}  // namespace ar
