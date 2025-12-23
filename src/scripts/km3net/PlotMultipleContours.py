#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import chi2
from matplotlib.lines import Line2D

# =====================
# File + style settings
# =====================

colors = ['teal','steelblue','blue','darkblue','firebrick','red','darkorange','gold'] 
files_config = [
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/sensitivity_study/final2/efstat_grid_results_withzenithcut0.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "teal",
        "linestyle": "-",
        "alpha": 1.0,
        "label": "Zenith = 0°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/sensitivity_study/final2/efstat_grid_results_withzenithcut5.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "steelblue",
        "linestyle": "-",
        "alpha": 0.9,
        "label": "Zenith = 5°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/sensitivity_study/final2/efstat_grid_results_withzenithcut10.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "blue",
        "linestyle": "-",
        "alpha": 0.9,
        "label": "Zenith = 10°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/sensitivity_study/final2/efstat_grid_results_withzenithcut15.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "darkblue",
        "linestyle": "-",
        "alpha": 0.9,
        "label": "Zenith = 15°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/sensitivity_study/final2/efstat_grid_results_withzenithcut30.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "firebrick",
        "linestyle": "-",
        "alpha": 0.9,
        "label": "Zenith = 30°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/sensitivity_study/final2/efstat_grid_results_withzenithcut60.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "red",
        "linestyle": "-",
        "alpha": 0.9,
        "label": "Zenith = 60°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/sensitivity_study/final2/efstat_grid_results_withzenithcut90.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "darkorange",
        "linestyle": "-",
        "alpha": 0.9,
        "label": "Zenith = 90°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/sensitivity_study/final2/efstat_grid_results_withzenithcut120.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "gold",
        "linestyle": "-",
        "alpha": 0.9,
        "label": "Zenith = 120°"
    },
]
files_configv = [
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final2/efstat_grid_results_withzenithcut0.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "teal",
        "linestyle": "-",
        "alpha": 1.0,
        "detector" : "arca",
        "label": "Minimum Zenith Angle = 0°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final2/efstat_grid_results_withzenithcut90.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "steelblue",
        "linestyle": "-",
        "alpha": 1.0,
        "detector" : "arca",
        "label": "Minimum Zenith Angle = 90°"
    }]

files_configc = [
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final2/efstat_grid_results_withzenithcut0.0_classicflux_Crab_30.0deg_1iterations.csv",
        "color": "teal",
        "linestyle": "-",
        "alpha": 1.0,
        "detector" : "arca",
        "label": "Minimum Zenith Angle = 0°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final2/efstat_grid_results_withzenithcut90.0_classicflux_Crab_30.0deg_1iterations.csv",
        "color": "steelblue",
        "linestyle": "-",
        "alpha": 1.0,
        "detector" : "arca",
        "label": "Minimum Zenith Angle = 90°"
    }]

files_config2 = [
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final4/efstat_grid_results_withzenithcut0.0_boundedflux_Vela_5.0deg_1iterations.csv",
        "color": "teal",
        "linestyle": "-",
        "alpha": 0.9,
        "detector" : "arca",
        "label": "Minimum Zenith Angle = 0"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final4/efstat_grid_results_withzenithcut30.0_boundedflux_Vela_5.0deg_1iterations.csv",
        "color": "steelblue",
        "linestyle": "-",
        "alpha": 0.9,
        "detector" : "arca",
        "label": "Minimum Zenith Angle = 30°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final4/efstat_grid_results_withzenithcut60.0_boundedflux_Vela_5.0deg_1iterations.csv",
        "color": "blue",
        "linestyle": "-",
        "alpha": 1.0,
        "detector" : "arca",
        "label": "Minimum Zenith Angle = 60°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final4/efstat_grid_results_withzenithcut90.0_boundedflux_Vela_5.0deg_1iterations.csv",
        "color": "darkblue",
        "linestyle": "-",
        "alpha": 1.0,
        "detector" : "arca",
        "label": "Minimum Zenith Angle = 90°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final4/efstat_grid_results_withzenithcut120.0_boundedflux_Vela_5.0deg_1iterations.csv",
        "color": "firebrick",
        "linestyle": "-",
        "alpha": 1.0,
        "detector" : "arca",
        "label": "Minimum Zenith Angle = 120°"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final4/efstat_grid_results_withzenithcut150.0_boundedflux_Vela_5.0deg_1iterations.csv",
        "color": "red",
        "linestyle": "-",
        "alpha": 1.0,
        "detector" : "arca",
        "label": "Minimum Zenith Angle = 150°"
    }
    
    ]


files_config2 = [
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final2/efstat_grid_results_withzenithcut0.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "steelblue",
        "linestyle": "-",
        "alpha": 1,
        "detector" : "arca",
        "label" : "ARCA"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final2/efstat_grid_results_withzenithcut0.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "teal",
        "linestyle": "-",
        "alpha": 1,
        "detector" : "orca",
        "label" : "ORCA"
    },
    {
        "path": "/home/hpc/capn/capn107h/software/psr/src/scripts/sensitivity_study/final2/efstat_grid_results_withzenithcut0.0_classicflux_Vela_30.0deg_1iterations.csv",
        "color": "darkorange",
        "linestyle": "--",
        "alpha": 1,
        "detector" : "combined",
        "label" : "Combined"
    }]

mpl.rcParams['axes.labelsize'] = 12
#mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['legend.fontsize'] = 11
mpl.rcParams['figure.titlesize'] = 14


# Output file
output_file = "contours_overlay_zenithcuts_Crab_30deg_classic.png"

# Epoch folding bins
nbins = 32

# Detector to plot
#detector_name = "arca"   # <--- change if you want another one


def plot_contours(files_config, nbins, output_file):
    dof = nbins - 1
    p_value = 2.87e-7  # 5 sigma
    chi2_threshold = chi2.ppf(1 - p_value, df=dof)
    print("chi2_threshold:", chi2_threshold)

    fig, ax = plt.subplots(figsize=(7, 6))
    legend_lines = []

    for color, cfg in zip(colors,files_config):
        f = cfg["path"]
        df = pd.read_csv(f)
        det = cfg["detector"]

        if "detector" not in df.columns:
            print(f"⚠️ Skipping {f}: no 'detector' column")
            continue

        df_sub = df[df['detector'].str.lower() == det]
        if df_sub.empty:
            print(f"⚠️ Skipping {f}: no data for detector '{det}'")
            continue

        pivot = df_sub.pivot(index='phi0', columns='gamma', values='mean_max_efstat')
        if pivot.empty:
            print(f"⚠️ Skipping {f}: pivot empty for detector '{det}'")
            continue

        data = pivot.to_numpy()
        gamma_grid = pivot.columns.to_numpy()
        phi0_grid = pivot.index.to_numpy()
        X, Y = np.meshgrid(gamma_grid, np.log10(phi0_grid))

        cs = ax.contour(
            X, Y, data,
            levels=[chi2_threshold],
            colors=[cfg["color"]],
            #colors=color,
            linestyles=cfg["linestyle"],
            alpha=cfg["alpha"],
            linewidths=2
        )

        legend_lines.append(Line2D([0], [0],
                                   color=cfg["color"],
                                   linestyle=cfg["linestyle"],
                                   alpha=cfg["alpha"],
                                   lw=2,
                                   label=cfg["label"]))

    # Axis formatting
    ax.set_xlabel(r"Spectral index $\gamma$")
    ax.set_ylabel(r"log10( Flux normalization $\Phi_0$ ) [GeV$^{-1}$ m$^{-2}$ s$^{-1}$]")
    ax.grid(which='both', alpha=0.5)
    ax.legend(handles=legend_lines, loc='upper left')
    plt.title(f"5σ Contours for different Detectors")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Saved contour overlay plot to {output_file}")


if __name__ == "__main__":
    plot_contours(files_configc, nbins, output_file)
