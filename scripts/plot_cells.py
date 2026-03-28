import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import numpy as np


def plot_cells_with_material(csv_file):
    # Load CSV
    df = pd.read_csv(csv_file)

    # Required columns
    required_cols = {"x_center", "y_center", "width", "height", "material_id"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")

    # Get unique materials
    materials = sorted(df["material_id"].unique())
    n_materials = len(materials)

    # Create a colormap
    cmap = cm.get_cmap("Set2", n_materials)

    # Map material_id -> color
    material_to_color = {mat: cmap(i) for i, mat in enumerate(materials)}

    # Create plot
    fig, ax = plt.subplots()

    # Plot each cell
    for _, row in df.iterrows():
        x_c = row["x_center"]
        y_c = row["y_center"]
        w = row["width"]
        h = row["height"]
        mat = row["material_id"]

        # Bottom-left corner
        x_bl = x_c - w / 2
        y_bl = y_c - h / 2

        rect = Rectangle(
            (x_bl, y_bl),
            w,
            h,
            facecolor=material_to_color[mat],
            edgecolor="black",
            linewidth=1.0,
        )

        ax.add_patch(rect)

    # Aspect ratio
    ax.set_aspect("equal", adjustable="box")
    ax.autoscale_view()

    # Labels
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title("Cell Grid with Material IDs")

    # Create legend
    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor=material_to_color[mat]) for mat in materials
    ]
    legend_labels = [f"Material {mat}" for mat in materials]

    ax.legend(legend_handles, legend_labels, title="Materials", loc="upper right")

    plt.show()


if __name__ == "__main__":
    plot_cells_with_material("cell_geometry.csv")
