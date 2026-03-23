#ifndef ACKROYD_RIYAIT_SIMULATION_H_
#define ACKROYD_RIYAIT_SIMULATION_H_

#include <vector>

#include "rectangular_region.h"

namespace ar {
class Simulation {
 public:
  /**
   * @brief Construct a new Simulation object
   *
   * @param regions Vector of RectangularRegion objects
   * @param si_tolerance Source iteration tolerance
   */
  Simulation(const std::vector<RectangularRegion>& regions,
             const double si_tolerance = 1e-8);

 private:
  const std::vector<RectangularRegion>&
      regions_;                /// Vector of regions (passed by reference)
  const double si_tolerance_;  /// scattering iteration convergence tolerance

  std::vector<std::vector<Cell>> cells_;  /// 2D vector of cells

  /**
   * @brief Pulls all cells from all regions and places them into a 1D vector of
   * Cells
   *
   * @return std::vector<Cell>
   */
  std::vector<Cell> PullCellsFromRegions();

  /**
   * @brief Takes flattened vector of cells, then sorts and reshapes them s.t.
   * cells_[0][0] will index to the bottom left corner of domain and
   * cells[nx-1][ny-1] will index to the top right.
   *
   * @param flattened_cells
   * @return std::vector<Cell>
   */
  std::vector<std::vector<Cell>> SortCells(
      const std::vector<Cell>& flattened_cells);

  /**
   * @brief Exports results to a CSV output file for postprocessing
   *
   */
  void ExportToCSV();
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_SIMULATION_H_