#ifndef ACKROYD_RIYAIT_SIMULATION_H_
#define ACKROYD_RIYAIT_SIMULATION_H_

#include <string>
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

  /**
   * @brief Exports all cells to a CSV file to be plotted and verify that the
   * geometry is correct before proceeding
   *
   * @param output_name Output file name (defaults to "cell_geometry.csv")
   */
  void ExportCellsToCSV(const std::string output_name = "cell_geometry.csv");

 private:
  const std::vector<RectangularRegion>&
      regions_;                /// Vector of regions (passed by reference)
  const double si_tolerance_;  /// scattering iteration convergence tolerance

  std::vector<std::vector<Cell>> cells_;  /// 2D vector of cells

  /**
   * @brief Pulls all cells from all regions and places them into a 1D vector
   * of Cells.
   *
   * @return std::vector<Cell>
   */
  std::vector<Cell> PullCellsFromRegions();

  /**
   * @brief Takes flattened vector of cells, then sorts and reshapes them s.t.
   * cells_[0][0] will index to the bottom left corner of the domain and
   * cells[nx-1][ny-1] will index to the top right.
   *
   * @param flattened_cells
   * @return std::vector<Cell>
   */
  std::vector<std::vector<Cell>> SortCells(
      const std::vector<Cell>& flattened_cells);

  int n_columns_;  /// Number of columns in 2D vector of cells
  int n_rows_;     /// Number of rows in 2D vector of cells

  /**
   * @brief Carries out transport sweep from Southwest corner to Northeast.
   * Returns a 2D vector of SCALAR flux contributions by multiplying the
   * cell-centered angular flux by the quadrature weight.
   *
   * @param quadrature_weight
   * @param x_cosine
   * @param y_cosine
   * @return std::vector<std::vector<double>>
   */
  std::vector<std::vector<double>> SweepNorthEast(
      const double quadrature_weight, const double x_cosine,
      const double y_cosine);

  /**
   * @brief Carries out transport sweep from Southeast corner to Northwest.
   * Returns a 2D vector of SCALAR flux contributions by multiplying the
   * cell-centered angular flux by the quadrature weight.
   *
   * @param quadrature_weight
   * @param x_cosine
   * @param y_cosine
   * @return std::vector<std::vector<double>>
   */
  std::vector<std::vector<double>> SweepNorthWest(
      const double quadrature_weight, const double x_cosine,
      const double y_cosine);

  /**
   * @brief Carries out transport sweep from Northwest corner to Southeast.
   * Returns a 2D vector of SCALAR flux contributions by multiplying the
   * cell-centered angular flux by the quadrature weight.
   *
   * @param quadrature_weight
   * @param x_cosine
   * @param y_cosine
   * @return std::vector<std::vector<double>>
   */
  std::vector<std::vector<double>> SweepSouthEast(
      const double quadrature_weight, const double x_cosine,
      const double y_cosine);

  /**
   * @brief Carries out transport sweep from Northeast corner to Southwest.
   * Returns a 2D vector of SCALAR flux contributions by multiplying the
   * cell-centered angular flux by the quadrature weight.
   *
   * @param quadrature_weight
   * @param x_cosine
   * @param y_cosine
   * @return std::vector<std::vector<double>>
   */
  std::vector<std::vector<double>> SweepSouthWest(
      const double quadrature_weight, const double x_cosine,
      const double y_cosine);

  /**
   * @brief Evaluates the cell-centered angular flux for a single cell.
   *
   * @param x_cosine
   * @param y_cosine
   * @param cell
   * @return double
   */
  double SweepStep(const double x_cosine, const double y_cosine, Cell& cell);
  /**
   * @brief Exports results to a CSV output file for postprocessing.
   *
   */
  void ExportResultsToCSV();
};
}  // namespace ar

#endif  // ACKROYD_RIYAIT_SIMULATION_H_
