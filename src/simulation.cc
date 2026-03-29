#include "simulation.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace ar {
Simulation::Simulation(const std::vector<RectangularRegion>& regions,
                       const double si_tolerance)
    : regions_(regions), si_tolerance_(si_tolerance) {
  std::vector<Cell> flattened_cells = PullCellsFromRegions();
  cells_ = SortCells(flattened_cells);
  n_columns_ = cells_[0].size();
  n_rows_ = cells_.size();
  ExportCellsToCSV();
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
  // Sort all indices by X
  std::vector<size_t> indices(flattened_cells.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
            [&flattened_cells](size_t a, size_t b) {
              return flattened_cells[a].cell_center().X() <
                     flattened_cells[b].cell_center().X();
            });

  std::vector<std::vector<Cell>> cells;
  std::vector<size_t> temp_indices;
  double curr_x = flattened_cells[indices[0]].cell_center().X();

  auto flush_group = [&]() {
    // Sort temp_indices by Y, then copy-construct Cells into a new row
    std::sort(temp_indices.begin(), temp_indices.end(),
              [&flattened_cells](size_t a, size_t b) {
                return flattened_cells[a].cell_center().Y() <
                       flattened_cells[b].cell_center().Y();
              });
    std::vector<Cell> row;
    for (size_t i : temp_indices) {
      row.push_back(flattened_cells[i]);  // copy-constructs
    }
    cells.push_back(std::move(row));
    temp_indices.clear();
  };

  for (size_t idx : indices) {
    double x = flattened_cells[idx].cell_center().X();
    if (x != curr_x) {
      flush_group();
      curr_x = x;
    }
    temp_indices.push_back(idx);
  }
  flush_group();  // handle the last group

  return cells;
}

void Simulation::ExportCellsToCSV(const std::string output_name) {
  std::ofstream output_file(output_name);

  // Check if the file was opened successfully
  if (!output_file.is_open()) {
    throw std::runtime_error(
        "Output file could not be opened in ExportCellsToCSV!");
  }

  output_file << "x_center,y_center,width,height,material_id" << "\n";

  for (auto& row : cells_) {
    for (auto& cell : row)
      output_file << cell.cell_center().X() << "," << cell.cell_center().Y()
                  << "," << cell.dx() << "," << cell.dy() << ","
                  << cell.material().id() << "\n";
  }
  output_file.close();
}

void Simulation::ExportResultsToCSV(const std::string output_name) {
  std::ofstream output_file(output_name);

  // Check if the file was opened successfully
  if (!output_file.is_open()) {
    throw std::runtime_error(
        "Output file could not be opened in ExportResultstoCSV!");
  }

  output_file << "x_center,y_center,center_flux" << "\n";

  for (auto& row : cells_) {
    for (auto& cell : row)
      output_file << cell.cell_center().X() << "," << cell.cell_center().Y()
                  << "," << cell.scalar_flux() << "\n";
  }
  output_file.close();
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
