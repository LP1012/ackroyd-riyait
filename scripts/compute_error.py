import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler

plt.style.use("matplotlib.rc")
plt.rcParams["axes.prop_cycle"] = cycler(
    color=mpl.colormaps["cividis"](np.linspace(0, 1, 11, endpoint=True))
)


def ComputeL2Error(true_solution_csv, test_csv):
    df_known = pd.read_csv(true_solution_csv)
    df_test = pd.read_csv(test_csv)

    df_known = df_known.sort_values(["x_center", "y_center"])
    df_test = df_test.sort_values(["x_center", "y_center", "volume"])

    xs_known = np.array(sorted(df_known["x_center"].unique()))
    ys_known = np.array(sorted(df_known["y_center"].unique()))

    xs_test = np.array(sorted(df_test["x_center"].unique()))
    ys_test = np.array(sorted(df_test["y_center"].unique()))
    volumes_test = df_test.pivot(
        index="y_center", columns="x_center", values="volume"
    ).to_numpy()

    flux_known = df_known.pivot(
        index="y_center", columns="x_center", values="center_flux"
    ).to_numpy()
    flux_test = df_test.pivot(
        index="y_center", columns="x_center", values="center_flux"
    )

    interpolator = RegularGridInterpolator(
        (ys_known, xs_known), flux_known, method="linear"
    )

    # Create all (y, x) combinations
    yy, xx = np.meshgrid(ys_test, xs_test, indexing="ij")
    points = np.column_stack([yy.ravel(), xx.ravel()])

    # Evaluate interpolator at all test points
    flux_interpolated = interpolator(points).reshape(len(ys_test), len(xs_test))

    error = flux_interpolated - flux_test.to_numpy()
    error_flattened = error.ravel()

    l2_error = np.sqrt(np.sum(error**2 * volumes_test))
    cell_volume = volumes_test[0, 0]
    print(f"L2-Error = {l2_error}, cell volume = {cell_volume}")

    return l2_error, cell_volume


def Plot(volumes, l2_errors):
    # create reference curves
    endpoints = np.array([volumes[0], volumes[-1]])
    plt.figure()
    plt.plot(
        endpoints,
        endpoints**2 * (l2_errors[0] / volumes[0] ** 2),
        label=r"$\mathcal{O}\left(V^2\right)$",
        linestyle="dashed",
    )
    plt.plot(
        endpoints,
        endpoints * (l2_errors[0] / volumes[0]),
        label=r"$\mathcal{O}\left(V\right)$",
        linestyle="dashed",
    )

    plt.plot(volumes, l2_errors, marker="o", label="Data")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Cell Volume")
    plt.ylabel(r"Discrete $L^2$ Norm of Error")
    plt.title(r"Grid Convergence Study for $N=32^2$ Ordinates")
    plt.legend()
    plt.savefig("grid_convergence.png")
    plt.close()


if __name__ == "__main__":
    l2_error_1, volume_1 = ComputeL2Error(
        "cell_fluxes_test_16.0.csv", "cell_fluxes_test_1.0.csv"
    )
    l2_error_2, volume_2 = ComputeL2Error(
        "cell_fluxes_test_16.0.csv", "cell_fluxes_test_2.0.csv"
    )
    l2_error_4, volume_4 = ComputeL2Error(
        "cell_fluxes_test_16.0.csv", "cell_fluxes_test_4.0.csv"
    )
    l2_error_8, volume_8 = ComputeL2Error(
        "cell_fluxes_test_16.0.csv", "cell_fluxes_test_8.0.csv"
    )

    errors = np.array([l2_error_1, l2_error_2, l2_error_4, l2_error_8])
    volumes = np.array([volume_1, volume_2, volume_4, volume_8])
    Plot(volumes, errors)
