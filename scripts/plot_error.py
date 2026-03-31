import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import linregress


def plot_error(csv_file):
    df = pd.read_csv(csv_file)

    required_cols = {"n_cells", "l2_norm"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")
    
    fig, ax = plt.subplots()

    xs = np.array(df["n_cells"])
    ys = np.zeros(len(df))

    last_l2norm = 1
    for i, row in df.iterrows():
        err = np.abs(row["l2_norm"] - last_l2norm)/row["l2_norm"]
        last_l2norm = row["l2_norm"]
        ys[i] = err

    ax.plot(xs, ys, marker='o')

    slope = linregress([np.log(x) for x in xs[1:]],
                       [np.log(y) for y in ys[1:]])
    
    slope_xs = np.linspace(256, 1e5, 500)
    ax.plot(slope_xs, [x**slope[0]*np.exp(slope[1]) for x in slope_xs],
            label=f'slope = {round(slope[0], 3)}')
    print(f"m = {slope[0]}, b = {slope[1]}")
    
    ax.set_xlabel("Number of Cells")
    ax.set_ylabel("Relative Error")
    ax.loglog()
    ax.grid()
    ax.legend()

    plt.tight_layout()
    plt.savefig("l2_norm.png", dpi=300)

if __name__ == "__main__":
    plot_error("l2_norm.csv")
