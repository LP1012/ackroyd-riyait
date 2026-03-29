import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def plot_cell_fluxes(csv_file):
    # Load CSV
    df = pd.read_csv(csv_file)

    # Required columns
    required_cols = {"x_center", "y_center", "center_flux"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")
    
    # Create plot
    fig, ax = plt.subplots()

    # Rebuild cell mesh with fluxes
    xs = set()
    ys = set()
    flux_mesh = []
    curr_x = df.iloc[0]["x_center"]

    temp_row = []
    for _, row in df.iterrows():
        if row["x_center"] == curr_x:
            temp_row.append(row["center_flux"])
            xs.add(row["x_center"])
            ys.add(row["y_center"])
        else:
            curr_x = row["x_center"]
            flux_mesh.append(temp_row)
            temp_row = []

            temp_row.append(row["center_flux"])
            xs.add(row["x_center"])
            ys.add(row["y_center"])

    flux_mesh.append(temp_row)
    X, Y = np.meshgrid(np.array(list(xs)), np.array(list(ys)))

    # Plot mesh
    ax.pcolormesh(X, Y, flux_mesh, shading='gouraud', cmap='magma')

    # Labels
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title("Scalar Flux over Problem Domain")

    plt.show()


if __name__ == "__main__":
    plot_cell_fluxes("cell_fluxes.csv")
