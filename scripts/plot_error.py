import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


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

    ax.plot(xs, ys)
    
    ax.set_xlabel("Number of Cells")
    ax.set_ylabel("Relative Error")
    ax.loglog()
    ax.grid()

    plt.tight_layout()
    plt.savefig("l2_norm.png", dpi=300)

if __name__ == "__main__":
    plot_error("l2_norm.csv")
