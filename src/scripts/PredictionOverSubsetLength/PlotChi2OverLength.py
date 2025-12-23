import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy.optimize import curve_fit


def power_law(x, alpha, C):
    return C * x**alpha


def plot_mean_max_efstat(csv_file, output_file=None):
    # Load data
    df = pd.read_csv(csv_file)

    # Extract x and y
    x = df["length_days"].values
    y = df["mean_max_efstat"].values

    # Fit log-log linear model
    popt, pcov = curve_fit(power_law, x, y)
    alpha, C = popt
    print(f"Fitted relation: chi2 ≈ {C:.3e} * T^{alpha:.3f}")

    # Generate smooth curve for fit
    x_fit = np.logspace(np.log10(min(x)), np.log10(max(x)), 200)
    y_fit = power_law(x_fit, alpha, C)

    # Create plot
    plt.figure(figsize=(7, 5))
    plt.plot(x, y, marker="o", linestyle="-", color="darkblue", label="Epoch Folding Results")
    plt.plot(x_fit, y_fit,color="darkorange", linestyle="--", label=f"Fit: χ² ≈ {C:.1e}·T^{alpha:.2f}")

    # Labels and title
    plt.xlabel("Length of Dataset (days)", fontsize=12)
    plt.ylabel(r"Mean Maximum $\chi^2$", fontsize=12)
    plt.title("Epoch Folding Statistics over Length (Vela)", fontsize=14)
    #plt.yscale("log")
    #plt.xscale("log")

    # Grid, legend
    plt.grid(True, linestyle="--", which="both", alpha=0.6)
    plt.legend()

    # Save or show
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Plot saved as {output_file}")
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot mean_max_efstat vs length from CSV.")
    parser.add_argument("csv_file", help="Path to the CSV file")
    parser.add_argument("--output", "-o", help="Optional output image file (e.g. plot.png)")
    args = parser.parse_args()

    plot_mean_max_efstat(args.csv_file, args.output)
