import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors


def plot_cell_fluxes(csv_file):
    df = pd.read_csv(csv_file)

    required_cols = {"x_center", "y_center", "center_flux"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")

    # Sort to guarantee consistent grid ordering
    df = df.sort_values(["x_center", "y_center"])

    xs = np.array(sorted(df["x_center"].unique()))
    ys = np.array(sorted(df["y_center"].unique()))

    # Pivot into a 2D array: rows=y, cols=x (what pcolormesh expects)
    flux_mesh = df.pivot(
        index="y_center", columns="x_center", values="center_flux"
    ).values

    X, Y = np.meshgrid(xs, ys)

    fig, ax = plt.subplots()

    mesh = ax.imshow(
        flux_mesh,
        norm=mcolors.LogNorm(),
        cmap="magma",
        interpolation="bicubic",
        origin="lower",
        extent=[xs.min(), xs.max(), ys.min(), ys.max()],
        aspect="auto",
    )

    # Add colorbar — this also confirms the value range at a glance
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Scalar Flux")

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    # ax.set_title("Scalar Flux over Problem Domain")

    # ax.set_xlim(-10, 10)
    # ax.set_ylim(-1, 10)
    # # ax.set_xticks(xs)
    # # ax.set_yticks(ys)
    # ax.grid(color="white", linewidth=0.5, linestyle="--", alpha=0.5)

    ax.axis("off")
    cbar.remove()

    plt.tight_layout()
    plt.savefig("cell_fluxes.png", dpi=300, bbox_inches="tight", pad_inches=0)

    fig, ax = plt.subplots()

    diag_dist = np.diag(X)
    diag_flux = np.diag(flux_mesh)
    mid = int(len(diag_dist) / 2)

    ax.plot(diag_dist[mid - 1 :], diag_flux[mid - 1 :])
    ax.set_yscale("log")

    ax.set_xlabel("X-Y Diagonal")
    ax.set_ylabel("Scalar Flux")

    plt.tight_layout()
    plt.savefig("diagonal_flux.png", dpi=300)


if __name__ == "__main__":
    plot_cell_fluxes("cell_fluxes.csv")
