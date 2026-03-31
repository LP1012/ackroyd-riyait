import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd


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

    print(f"L2-Error = {l2_error}, cell volume = {volumes_test[0,0]}")


if __name__ == "__main__":
    ComputeL2Error("cell_fluxes_known.csv", "cell_fluxes_test.csv")
