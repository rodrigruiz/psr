import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
# =======================
# CONFIGURATION
# =======================
csv_files = [
    #"2025-08-23_17-10_efstat_grid_results_withzenithcut0.0_boundedflux_Vela_30.0deg.csv",
    #"2025-08-23_17-10_efstat_grid_results_withzenithcut5.0_boundedflux_Vela_30.0deg.csv",
    #"2025-08-23_17-10_efstat_grid_results_withzenithcut10.0_boundedflux_Vela_30.0deg.csv",
    #"2025-08-23_17-10_efstat_grid_results_withzenithcut15.0_boundedflux_Vela_30.0deg.csv",
    #"2025-08-23_17-10_efstat_grid_results_withzenithcut30.0_boundedflux_Vela_30.0deg.csv",
    #"2025-08-23_17-10_efstat_grid_results_withzenithcut60.0_boundedflux_Vela_30.0deg.csv",
    #"2025-08-23_17-09_efstat_grid_results_withzenithcut90.0_boundedflux_Vela_30.0deg.csv",
    "2025-08-23_19-45_efstat_grid_results_withzenithcut100.0_boundedflux_Vela_30.0deg.csv",
    "2025-08-23_19-45_efstat_grid_results_withzenithcut110.0_boundedflux_Vela_30.0deg.csv",
    "2025-08-23_17-10_efstat_grid_results_withzenithcut120.0_boundedflux_Vela_30.0deg.csv",
    "2025-08-23_19-45_efstat_grid_results_withzenithcut130.0_boundedflux_Vela_30.0deg.csv",
    "2025-08-23_19-45_efstat_grid_results_withzenithcut140.0_boundedflux_Vela_30.0deg.csv",
    "2025-08-23_19-45_efstat_grid_results_withzenithcut150.0_boundedflux_Vela_30.0deg.csv",
    "efstat_grid_results_withzenithcut160.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut170.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut180.0_boundedflux_Vela_30.0deg_1iterations.csv",
]

csv_files = [
    "efstat_grid_results_withzenithcut0.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut10.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut20.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut30.0_boundedflux_Vela_30.0deg_1iterations.csv", 
    "efstat_grid_results_withzenithcut40.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut50.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut60.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut70.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut80.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut90.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut100.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut110.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut120.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut130.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut140.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut150.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut160.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut170.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "efstat_grid_results_withzenithcut180.0_boundedflux_Vela_30.0deg_1iterations.csv",
]

csv_files = [
    "final2/efstat_grid_results_withzenithcut90.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "final2/efstat_grid_results_withzenithcut0.0_boundedflux_Vela_30.0deg_1iterations.csv",
    "final2/efstat_grid_results_withzenithcut90.0_boundedflux_Crab_30.0deg_1iterations.csv",
    "final2/efstat_grid_results_withzenithcut0.0_boundedflux_Crab_30.0deg_1iterations.csv",
]

colors = ["red", "green", "blue"]  # One color per zenith_cut
colors = ['grey','teal','steelblue','blue','darkblue','firebrick','red','darkorange','gold'] 
detector_to_plot = "arca"
sigma_threshold = 88.0  # 5σ detection
# Compute 5σ chi² threshold for nbins - 1 degrees of freedom
nbins=32
dof = nbins - 1
p_value = 2.87e-7  # 5 sigma
chi2_threshold = chi2.ppf(1 - p_value, df=dof)
print("chi2_threshold: ", chi2_threshold)
sigma_threshold = chi2_threshold
# =======================
# FUNCTIONS
# =======================
def find_threshold_flux_log(subdf, sigma_level=sigma_threshold):
    """
    Interpolate threshold phi0 in log-space.
    Expects subdf with columns: phi0, mean_max_efstat.
    """
    subdf = subdf.sort_values("phi0")
    phi = subdf["phi0"].values
    ef = subdf["mean_max_efstat"].values
    
    mask = np.isfinite(ef)
    phi = phi[mask]
    ef = ef[mask]
    
    if np.all(ef < sigma_level):
        return np.nan
    
    log_phi = np.log10(phi)
    thr_logphi = np.interp(sigma_level, ef, log_phi)
    return 10**thr_logphi


def compute_sensitivity(df, sigma_level=sigma_threshold):
    """
    Compute threshold phi0 for each (gamma, detector, angle, zenith_cut).
    Returns a tidy DataFrame.
    """
    results = []
    grouped = df.groupby(["gamma", "detector", "angle", "zenith_cut"])
    
    for (gamma, detector, angle, zenith_cut), subdf in grouped:
        thr = find_threshold_flux_log(subdf, sigma_level)
        results.append({
            "gamma": gamma,
            "detector": detector,
            "angle": angle,
            "zenith_cut": zenith_cut,
            "threshold_phi0": thr
        })
    
    return pd.DataFrame(results)


# =======================
# LOAD DATA
# =======================
all_dfs = [pd.read_csv(f) for f in csv_files]
df_all = pd.concat(all_dfs, ignore_index=True)

# Ensure angle and zenith_cut exist
if "angle" not in df_all.columns:
    df_all["angle"] = np.nan
if "zenith_cut" not in df_all.columns:
    df_all["zenith_cut"] = np.nan

# =======================
# COMPUTE THRESHOLDS
# =======================
df_sens = compute_sensitivity(df_all, sigma_level=sigma_threshold)

# =======================
# PLOTTING MULTIPLE ZENITH CUTS
# =======================
plt.figure(figsize=(8,6))

# Identify unique zenith cuts
unique_zeniths = sorted(df_sens["zenith_cut"].dropna().unique())

for i, zenith in enumerate(unique_zeniths):
    sub = df_sens[(df_sens["detector"] == detector_to_plot) &
                  (df_sens["zenith_cut"] == zenith)]
    
    if len(sub) == 0:
        continue
    
    # Sort by gamma for smooth lines
    sub = sub.sort_values("gamma")
    
    plt.plot(sub["gamma"], sub["threshold_phi0"], 
             color=colors[i % len(colors)], marker='o',
             label=f"Zenith cut {zenith}°")

plt.yscale("log")
plt.xlabel("Spectral index γ")
plt.ylabel("Threshold flux φ₀ (5σ)")
plt.title(f"Sensitivity curves vs γ ({detector_to_plot.upper()})")
plt.grid(True, which='both', ls='--', lw=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("Test.png")
