#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from scipy.stats import chi2

mpl.rcParams['axes.labelsize'] = 12
#mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['legend.fontsize'] = 1

def plot_detector_heatmaps(df, nbins, output_file, zenith, angle, fluxtype):
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
    chi2_threshold = chi2.ppf(1 - p_value, df=dof) #+ 10
    print("chi2_threshold: ", chi2_threshold)

    fig, axs = plt.subplots(1, len(detectors), figsize=(6 * len(detectors), 5), constrained_layout=True)

    # Ensure axs is iterable even if only one detector
    if len(detectors) == 1:
        axs = [axs]

    for ax, detector in zip(axs, detectors):
        sub = df[df['detector'] == detector]

        # Pivot data: rows=phi0, cols=gamma
        pivot = sub.pivot(index='phi0', columns='gamma', values='mean_max_efstat')
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
        ax.set_title(f"Epoch Folding Sensitivity ({detector}, {fluxtype}, {zenith}, {angle})")
        ax.set_xlabel(r"Spectral index $\gamma$")
        ax.set_ylabel(r"log10( Flux normalization $\Phi_0$ ) [GeV${}^{-1}$ m${}^{-2}$ s${}^{-1}$]")
        #ax.set_xticks(gamma_values)
        #ax.set_yticks(np.log10(phi0_values))
        #ax.set_yticklabels([f"{v:.0e}" for v in phi0_values])
        #ax.set_yscale('log')
        ax.grid(alpha=0.5)

        fig.colorbar(im, ax=ax, label="Maximum chi2")

    #plt.suptitle(fr"Maximum Chi2 over Phi0 and gamma per Detector with $5\sigma$ Contour", fontsize=16)
    plt.savefig(output_file, dpi=300)
    print(f"Saved plot to {output_file}")

def plot_arca_vs_combined(df, nbins, output_file):
    """
    Plot combined detector heatmap with ARCA-only and COMBINED 5σ contours.
    """
    phi0_values = sorted(df['phi0'].unique())
    gamma_values = sorted(df['gamma'].unique())

    dof = nbins - 1
    p_value = 2.87e-7  # 5 sigma
    chi2_threshold = chi2.ppf(1 - p_value, df=dof)
    print("chi2_threshold: ", chi2_threshold)
    fig, ax = plt.subplots(figsize=(7, 6))

    # --- Heatmap for COMBINED detector
    df_comb = df[df['detector'].str.lower() == "combined"]
    pivot_combined = df_comb.pivot(index='phi0', columns='gamma', values='mean_max_efstat')
    data_combined = pivot_combined.to_numpy()
    gamma_grid = pivot_combined.columns.to_numpy()
    phi0_grid = pivot_combined.index.to_numpy()

    im = ax.imshow(
        data_combined,
        origin='lower',
        aspect='auto',
        cmap='viridis',
        norm=LogNorm(vmin=np.nanmin(data_combined[data_combined > 0]), vmax=np.nanmax(data_combined)),
        extent=[min(gamma_values), max(gamma_values), np.log10(min(phi0_values)), np.log10(max(phi0_values))]
    )

    X, Y = np.meshgrid(gamma_grid, np.log10(phi0_grid))

    # --- 5σ contour for COMBINED
    cs_comb = ax.contour(X, Y, data_combined, levels=[chi2_threshold], colors='white', linewidths=2)
    cs_comb.collections[0].set_label("Combined 5σ")

    # --- 5σ contour for ARCA
    df_arca = df[df['detector'].str.lower() == "arca"]
    pivot_arca = df_arca.pivot(index='phi0', columns='gamma', values='mean_max_efstat')
    data_arca = pivot_arca.to_numpy()
    gamma_grid_arca = pivot_arca.columns.to_numpy()
    phi0_grid_arca = pivot_arca.index.to_numpy()
    X_arca, Y_arca = np.meshgrid(gamma_grid_arca, np.log10(phi0_grid_arca))

    cs_arca = ax.contour(X_arca, Y_arca, data_arca, levels=[chi2_threshold],
                         colors='grey', linewidths=2, linestyles="--")
    cs_arca.collections[0].set_label("ARCA 5σ")

    # --- Formatting
    ax.set_xlabel(r"Spectral index $\gamma$")
    ax.set_ylabel(r"Flux normalization $\Phi_0$ [GeV${}^{-1}$ m${}^{-2}$ s${}^{-1}$]")
    ax.set_xticks(gamma_values)
    ax.set_yticks(np.log10(phi0_values))
    ax.set_yticklabels([f"{v:.0e}" for v in phi0_values])
    fig.colorbar(im, ax=ax, label="Maximum chi2")
    from matplotlib.lines import Line2D

    # Create legend manually with proxy lines
    legend_lines = [
        Line2D([0], [0], color='white', lw=2, label='Combined 5σ'),
        Line2D([0], [0], color='grey', lw=2, linestyle='--', label='ARCA 5σ')
    ]
    ax.legend(handles=legend_lines, loc='upper right')


    plt.title("Comparison: ARCA vs COMBINED 5σ Sensitivity", fontsize=14)
    plt.savefig(output_file, dpi=300)
    print(f"Saved ARCA vs COMBINED plot to {output_file}")



def main():
    parser = argparse.ArgumentParser(description="Plot heatmaps of max_efstat with 5σ contour from CSV.")
    parser.add_argument("--csv", help="Path to the CSV file containing data")
    parser.add_argument("--nbins", type=int, default=32, help="Number of bins in epoch folding (default: 32)")
    parser.add_argument("--output", default="heatmaps2.png", help="Output PNG file name (default: heatmaps.png)")
    parser.add_argument("--zenith",type=float,default=90.0)
    parser.add_argument("--angle",type=float, default=30.0)
    parser.add_argument("--fluxtype", default="classic")
    args = parser.parse_args()

    # Load CSV
    df = pd.read_csv(args.csv)

    # Validate columns
    required_cols = {"detector", "phi0", "gamma", "mean_max_efstat","mean_max_efstat_error"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")

    # Create plots
    plot_detector_heatmaps(df, args.nbins, args.output, args.zenith, args.angle, args.fluxtype)
    #plot_arca_vs_combined(df, args.nbins, "heatmap_arca_vs_combined.png")


if __name__ == "__main__":
    main()
