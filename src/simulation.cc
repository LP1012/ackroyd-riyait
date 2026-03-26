#include "simulation.h"

#include <cmath>

namespace ar {
Simulation::Simulation(const std::vector<RectangularRegion>& regions,
                       const double si_tolerance)
    : regions_(regions), si_tolerance_(si_tolerance) {
  std::vector<Cell> flattened_cells = PullCellsFromRegions();
  cells_ = SortCells(flattened_cells);
  n_columns_ = cells_[0].size();
  n_rows_ = cells_.size();
}

std::vector<Cell> Simulation::PullCellsFromRegions() {
  std::vector<Cell> flattened_cells;
  for (auto& region : regions_) {
    for (auto& cell : region.cells()) flattened_cells.push_back(cell);
  }
  return flattened_cells;
}

std::vector<std::vector<Cell>> Simulation::SortCells(
    const std::vector<Cell>& flattened_cells) {
  std::vector<std::vector<Cell>> cells_;

  for (auto j = 0; j < n_columns_; j++) {
    std::vector<Cell> temp_;
    for (auto i = 0; i < n_rows_; i++) {
      temp_.push_back(flattened_cells[n_columns_ * i + j]);
    }
    cells_.push_back(temp_);
  }
}

std::vector<std::vector<double>> Simulation::SweepNorthEast(
    const double quadrature_weight, const double x_cosine,
    const double y_cosine) {
  std::vector<std::vector<double>> scalar_flux_contribution(
      n_rows_,
      std::vector<double>(n_columns_, 0.0));  // preallocates to all zeros

  for (auto j = 0; j < n_rows_; j++) {
    for (auto i = 0; i < n_columns_; i++) {
      double cell_center_flux = SweepStep(x_cosine, y_cosine, cells_[i][j]);
      scalar_flux_contribution[i][j] += cell_center_flux * quadrature_weight;
      cells_[i][j].SetCenterFlux(cell_center_flux);
      cells_[i][j].SetEastFlux();
      cells_[i][j].SetNorthFlux();
    }
  }
  return scalar_flux_contribution;
}

std::vector<std::vector<double>> Simulation::SweepNorthWest(
    const double quadrature_weight, const double x_cosine,
    const double y_cosine) {
  std::vector<std::vector<double>> scalar_flux_contribution(
      n_rows_,
      std::vector<double>(n_columns_, 0.0));  // preallocates to all zeros

  for (auto j = 0; j < n_rows_; j++) {
    for (auto i = n_columns_ - 1; i > -1; i--) {
      double cell_center_flux = SweepStep(x_cosine, y_cosine, cells_[i][j]);
      scalar_flux_contribution[i][j] += cell_center_flux * quadrature_weight;
      cells_[i][j].SetCenterFlux(cell_center_flux);
      cells_[i][j].SetWestFlux();
      cells_[i][j].SetNorthFlux();
    }
  }
  return scalar_flux_contribution;
}

std::vector<std::vector<double>> Simulation::SweepSouthEast(
    const double quadrature_weight, const double x_cosine,
    const double y_cosine) {
  std::vector<std::vector<double>> scalar_flux_contribution(
      n_rows_,
      std::vector<double>(n_columns_, 0.0));  // preallocates to all zeros

  for (auto j = n_rows_ - 1; j > -1; j--) {
    for (auto i = 0; i > n_columns_; i++) {
      double cell_center_flux = SweepStep(x_cosine, y_cosine, cells_[i][j]);
      scalar_flux_contribution[i][j] += cell_center_flux * quadrature_weight;
      cells_[i][j].SetCenterFlux(cell_center_flux);
      cells_[i][j].SetEastFlux();
      cells_[i][j].SetSouthFlux();
    }
  }
  return scalar_flux_contribution;
}

std::vector<std::vector<double>> Simulation::SweepSouthWest(
    const double quadrature_weight, const double x_cosine,
    const double y_cosine) {
  std::vector<std::vector<double>> scalar_flux_contribution(
      n_rows_,
      std::vector<double>(n_columns_, 0.0));  // preallocates to all zeros

  for (auto j = n_rows_ - 1; j > -1; j--) {
    for (auto i = n_columns_ - 1; i > -1; i--) {
      double cell_center_flux = SweepStep(x_cosine, y_cosine, cells_[i][j]);
      scalar_flux_contribution[i][j] += cell_center_flux * quadrature_weight;
      cells_[i][j].SetCenterFlux(cell_center_flux);
      cells_[i][j].SetWestFlux();
      cells_[i][j].SetSouthFlux();
    }
  }
  return scalar_flux_contribution;
}

double Simulation::SweepStep(const double x_cosine, const double y_cosine,
                             Cell& cell) {
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
