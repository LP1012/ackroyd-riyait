#include "simulation.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace ar {
Simulation::Simulation(const std::vector<RectangularRegion>& regions,
                       const int n_ordinates, const double si_tolerance)
    : regions_(regions),
      spherical_quadrature_(GaussLegendreChebyshev(n_ordinates)),
      si_tolerance_(si_tolerance) {
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

double Simulation::ScalarFluxL2Norm() {
  double sum = 0;
  for (auto& cell_row : cells_) {
    for (auto& cell : cell_row) sum += cell.ScalarFluxL2();
  }
  return sum;
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

void Simulation::Run() {
  printf("Beginning simulation...\n");
  printf("  number of cells = %d\n", n_rows_ * n_columns_);
  printf("  scattering source iteration tolerance = %.3e\n", si_tolerance_);
  printf("  number of ordinates per quadrature set = %d\n\n",
         spherical_quadrature_.n_ordinates());

  printf("Beginning scattering iterations...\n");
  ScatteringIteration();
}
void Simulation::ScatteringIteration() {
  double relative_error = 1.0;  // starting dummy value

  double old_scalar_flux_l2 = ScalarFluxL2Norm();

  unsigned int count = 1;
  while (relative_error > si_tolerance_ && count <= 5) {
    InitializeCells();

    SweepCells();

    double new_scalar_flux_l2 = ScalarFluxL2Norm();
    relative_error =
        std::abs(new_scalar_flux_l2 - old_scalar_flux_l2) / new_scalar_flux_l2;
    printf("  Iteration : %4d, Relative Error = %.3e\n", count, relative_error);

    count++;
    old_scalar_flux_l2 = new_scalar_flux_l2;
  }
}

void Simulation::InitializeCells() {
  for (auto& cell_row : cells_) {
    for (auto& cell : cell_row) {
      cell.SetCellSource();
      cell.ClearScalarFlux();
    };
  }
}

void Simulation::SweepCells() {
  for (auto& triplet : spherical_quadrature_.GetTriples()) {
    // printf("    mu = %.4f, eta = %.4f\n", triplet.mu, triplet.eta);
    for (auto& cell_row : cells_) {
      for (auto& cell : cell_row) cell.ResetBoundaryFluxes();
    }
    if (triplet.mu > 0 && triplet.eta > 0)
      SweepNorthEast(triplet);
    else if (triplet.mu > 0 && triplet.eta < 0)
      SweepSouthEast(triplet);
    else if (triplet.mu < 0 && triplet.eta > 0)
      SweepNorthWest(triplet);
    else if (triplet.mu < 0 && triplet.eta < 0)
      SweepSouthWest(triplet);
    else
      throw std::runtime_error(
          "Direction of sweeping not able to be determined!");
  }
}

void Simulation::SweepNorthEast(const GLCTriplet& glc_triplet) {
  // printf("    Sweeping North East...");
  for (auto j = 0; j < n_rows_; j++) {
    for (auto i = 0; i < n_columns_; i++) {
      Cell& current_cell = cells_[i][j];

      current_cell.SetCenterFlux(
          SweepStep(glc_triplet.mu, glc_triplet.eta, current_cell));
      current_cell.AddPartialScalarFlux(current_cell.center_flux() *
                                        glc_triplet.weight);
      if (i + 1 < n_columns_)
        cells_[i + 1][j].SetWestFlux(2.0 * current_cell.center_flux() -
                                     current_cell.west_flux());
      if (j + 1 < n_rows_)
        cells_[i][j + 1].SetSouthFlux(2.0 * current_cell.center_flux() -
                                      current_cell.south_flux());
    }
  }
  // printf("Done.\n");
}

void Simulation::SweepNorthWest(const GLCTriplet& glc_triplet) {
  // printf("    Sweeping North West...");
  for (auto j = 0; j < n_rows_; j++) {
    for (auto i = n_columns_ - 1; i > -1; i--) {
      Cell& current_cell = cells_[i][j];

      current_cell.SetCenterFlux(
          SweepStep(glc_triplet.mu, glc_triplet.eta, current_cell));
      current_cell.AddPartialScalarFlux(current_cell.center_flux() *
                                        glc_triplet.weight);

      if (i - 1 >= 0)
        cells_[i - 1][j].SetEastFlux(2.0 * current_cell.center_flux() -
                                     current_cell.east_flux());
      if (j + 1 < n_rows_)
        cells_[i][j + 1].SetSouthFlux(2.0 * current_cell.center_flux() -
                                      current_cell.south_flux());
    }
  }
  // printf("Done.\n");
}

void Simulation::SweepSouthEast(const GLCTriplet& glc_triplet) {
  // printf("    Sweeping South East...");
  for (auto j = n_rows_ - 1; j > -1; j--) {
    for (auto i = 0; i > n_columns_; i++) {
      Cell& current_cell = cells_[i][j];

      current_cell.SetCenterFlux(
          SweepStep(glc_triplet.mu, glc_triplet.eta, current_cell));
      current_cell.AddPartialScalarFlux(current_cell.center_flux() *
                                        glc_triplet.weight);

      if (i + 1 < n_columns_)
        cells_[i + 1][j].SetWestFlux(2.0 * current_cell.center_flux() -
                                     current_cell.west_flux());
      if (j - 1 >= 0)
        cells_[i][j - 1].SetNorthFlux(2.0 * current_cell.center_flux() -
                                      current_cell.north_flux());
    }
  }
  // printf("Done.\n");
}

void Simulation::SweepSouthWest(const GLCTriplet& glc_triplet) {
  // printf("    Sweeping South West...");
  for (auto j = n_rows_ - 1; j > -1; j--) {
    for (auto i = n_columns_ - 1; i > -1; i--) {
      Cell& current_cell = cells_[i][j];

      current_cell.SetCenterFlux(
          SweepStep(glc_triplet.mu, glc_triplet.eta, current_cell));
      current_cell.AddPartialScalarFlux(current_cell.center_flux() *
                                        glc_triplet.weight);
      if (i - 1 >= 0)
        cells_[i - 1][j].SetEastFlux(2.0 * current_cell.center_flux() -
                                     current_cell.east_flux());
      if (j - 1 >= 0)
        cells_[i][j - 1].SetNorthFlux(2.0 * current_cell.center_flux() -
                                      current_cell.north_flux());
    }
  }
  // printf("Done.\n");
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

  double value = numerator / denominator;
  return (!std::isnan(value)) ? value : 0;
}

}  // namespace ar
