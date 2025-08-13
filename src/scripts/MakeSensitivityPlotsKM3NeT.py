#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import chi2

def plot_detector_heatmaps(df, nbins, output_file):
    """
    Creates three heatmaps of max_efstat (one per detector),
    with gamma on x-axis and phi0 on y-axis.
    Adds a contour line showing the 5σ detection threshold for given nbins.
    Saves the plot to a PNG file.
    """
    detectors = df['detector'].unique()
    phi0_values = sorted(df['phi0'].unique())
    gamma_values = sorted(df['gamma'].unique())

    # Compute 5σ chi² threshold for nbins - 1 degrees of freedom
    dof = nbins - 1
    p_value = 2.87e-7  # 5 sigma
    chi2_threshold = chi2.ppf(1 - p_value, df=dof)

    fig, axs = plt.subplots(1, len(detectors), figsize=(6 * len(detectors), 5), constrained_layout=True)

    # Ensure axs is iterable even if only one detector
    if len(detectors) == 1:
        axs = [axs]

    for ax, detector in zip(axs, detectors):
        sub = df[df['detector'] == detector]

        # Pivot data: rows=phi0, cols=gamma
        pivot = sub.pivot(index='phi0', columns='gamma', values='max_efstat')
        data = pivot.to_numpy()
        gamma_grid = pivot.columns.to_numpy()
        phi0_grid = pivot.index.to_numpy()

        # Mask invalid data
        data = np.where(pd.isna(data), np.nan, data)

        # Heatmap
        im = ax.imshow(
            data,
            origin='lower',
            aspect='auto',
            cmap='viridis',
            norm=LogNorm(vmin=np.nanmin(data[data > 0]), vmax=np.nanmax(data)),
            extent=[
                min(gamma_values), max(gamma_values),
                np.log10(min(phi0_values)), np.log10(max(phi0_values))
            ]
        )

        # Overlay 5σ contour
        X, Y = np.meshgrid(gamma_grid, np.log10(phi0_grid))
        cs = ax.contour(
            X, Y, data,
            levels=[chi2_threshold],
            colors='white',
            linewidths=1.5
        )
        ax.clabel(cs, fmt={chi2_threshold: "5σ"}, inline=True, fontsize=10)

        # Axis formatting
        ax.set_title(f"Detector: {detector}")
        ax.set_xlabel("γ (gamma)")
        ax.set_ylabel("log10(φ₀)")
        ax.set_xticks(gamma_values)
        ax.set_yticks(np.log10(phi0_values))
        ax.set_yticklabels([f"{v:.0e}" for v in phi0_values])

        fig.colorbar(im, ax=ax, label="max_efstat (log scale)")

    plt.suptitle(f"max_efstat Heatmaps per Detector with 5σ Contour (nbins={nbins})", fontsize=16)
    plt.savefig(output_file, dpi=300)
    print(f"Saved plot to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Plot heatmaps of max_efstat with 5σ contour from CSV.")
    parser.add_argument("--csv_file", help="Path to the CSV file containing data")
    parser.add_argument("--nbins", type=int, default=32, help="Number of bins in epoch folding (default: 32)")
    parser.add_argument("--output", default="heatmaps.png", help="Output PNG file name (default: heatmaps.png)")
    args = parser.parse_args()

    # Load CSV
    df = pd.read_csv(args.csv_file)

    # Validate columns
    required_cols = {"detector", "phi0", "gamma", "max_efstat"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")

    # Create plots
    plot_detector_heatmaps(df, args.nbins, args.output)

if __name__ == "__main__":
    main()
