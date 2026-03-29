import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os

cwd = os.getcwd()

def plot_cell_fluxes(csv_file):
    # Load CSV
    df = pd.read_csv(csv_file)

    # Required columns
    required_cols = {"x_center", "y_center", "center_flux"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")
    
    # Create plot
    fig, ax = plt.subplots()

    # Plot flux
    coords = np.array([])
    flux_mesh = np.array([])
    curr_x = df.iloc[0]["x_center"]

    temp_row = np.array([])
    for _, row in df.iterrows():
        if row["x_center"] == curr_x:
            np.append(temp_row)

    X, Y = np.meshgrid(xs, ys)

    ax.contourf(X, Y, flx)

    plt.show()


if __name__ == "__main__":
    plot_cell_fluxes("cell_fluxes.csv")
