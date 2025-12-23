#!/usr/bin/env python3
"""
plot_injection_diagnostics.py

Generate diagnostic plots for MVMD and sinusoidal signal injections.
- Plots the underlying counts function for chosen parameters
- Plots sampled events (vertical lines) vs. counts
- Saves results into a specified output directory

Usage:
    python plot_injection_diagnostics.py --output plots --pulseshape mvm --method classic \
        --frequency 2.0 --phi 0.0 --a 5.0 --baseline 1.0 --kappa 2.0 --num_periods 10
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import matplotlib as mpl
#mpl.rcParams['font.size'] = 11
mpl.rcParams['axes.labelsize'] = 12
#mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['legend.fontsize'] = 1
# --- Import your functions ---
#from your_module import MVMD, sinusoid   # <-- replace with actual module name


def MVMD(t, f, phi, kappa, a, baseline=None):
    """Modivied von Mises distribution (MVMD).
        Reference: "Fourier Techniques for Very Long Astrophysical Time-Series Analysis" 
                    Scott M. Ransom et al 2002 AJ 124 1788
                    
    Parameters
    ----------
        t : array_like
            Values to evaluate the MVMD at.
        kappa : float
            Shape parameter giving the width of the function.
        a : float
            Amplitude of the Pulse. Equates to the area of one pulse.
        f : float
            Frequency of the pulse train.
        phi : float
            Phase offset of the pulse train.
    
    Returns
    -------
        y : 1darray
            Evaluated values of the MVMD.
    
    For kappa -> 0 the MVMD converges towards a sinusod.
    For kappa -> infinity the MVMD converges towards a Gaussian. 1/kappa corresponds to sigma^2. Keep in mind, that  
        
    """
    
    y = a * ( np.exp(kappa*np.cos(2*np.pi*f*t+phi))-np.exp(-kappa))/(np.i0(kappa) - np.exp(-kappa))
    
    if baseline is not None:
        y += baseline
        
    return y

def sinusoid(times, frequency, baseline, amplitude, phase):
    return baseline + amplitude * np.sin(2 * np.pi * (frequency * times + phase))

def plot_signal_examples(
        output_dir,
        pulseshape,
        method,
        frequency,
        phi,
        a,
        baseline=0.0,
        kappa=None,
        bin_time=1.0,
        num_periods=5,
        high_res_dt=0.01,
        random_seed=42
    ):
    """
    Generate and save diagnostic plots for injected signals.
    """

    np.random.seed(random_seed)
    period = 1.0 / frequency
    t_max = num_periods * period

    # --- High resolution timeline ---
    times = np.arange(0, t_max, high_res_dt)

    # --- Generate underlying counts ---
    if pulseshape == 'mvm':
        counts = MVMD(times, frequency, phi, kappa, a, baseline=baseline)
    elif pulseshape == 'sine':
        counts = sinusoid(times, frequency, baseline, a, phi)
    else:
        raise ValueError("Invalid pulseshape")

    # --- Simulate events ---
    if method == 'classic':
        prob = counts / np.sum(counts)
        n_events = int(num_periods * a)  # heuristic number of events
        sampled_times = np.random.choice(times, size=n_events, p=prob)
    elif method == 'base':
        # Bin the curve
        binned_time = np.arange(0, t_max, bin_time)
        if pulseshape == 'mvm':
            binned_counts = MVMD(binned_time, frequency, phi, kappa, a, baseline=baseline)
        else:
            binned_counts = sinusoid(binned_time, frequency, baseline, a, phi)
        prob = binned_counts / np.sum(binned_counts)
        n_events = int(num_periods * a)
        sampled_times = np.random.choice(binned_time, size=n_events, p=prob)
    else:
        raise ValueError("Invalid method")


    fig, ax1 = plt.subplots(figsize=(7, 5))

    # Primary axis: underlying signal
    
    ax1.vlines(sampled_times, ymin=min(counts)*0.9, ymax=max(counts)*1.1,
               color="grey", alpha=0.6, lw=0.8, label="Sampled events")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Counts / underlying function")
    ax1.tick_params(axis="y")

    # Secondary axis: histogram of events
    ax2 = ax1.twinx()
    hist_bins = int(num_periods * frequency * 20)  # heuristic: 20 bins per period
    hist_vals, bin_edges = np.histogram(sampled_times, bins=hist_bins, range=(0, t_max))
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    ax2.step(bin_centers, hist_vals, where="mid", color="darkorange", alpha=0.7, label="Sampled event histogram")
    ax2.set_ylabel("Sampled events per bin")
    ax2.tick_params(axis="y")

    ax1.plot(times, counts, label="Counts function", color="darkblue")

    # Title and legend
    ax1.set_title(f"Event Time Sampling ({pulseshape}, {method}, f={frequency:.2e}Hz, kappa={kappa:.2f}, {n_events} events)")
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    #ax1.legend(lines1 + lines2, labels1 + labels2,framealpha=0.9, fontsize="small", loc="lower left")
    legend = fig.legend(
        lines1 + lines2, labels1 + labels2,
        loc="upper left",
        bbox_to_anchor=(0.02, 0.98),   # relative to axes area
        bbox_transform=ax1.transAxes,  # <-- use axes fraction coords
        fontsize="small", framealpha=0.9
        )

    fname_events = os.path.join(
        output_dir,
        f"events_hist_{pulseshape}_{method}_f{frequency}_T{t_max:.2f}_kappa{kappa}_{n_events}events_75fmt.png"
    )
    plt.savefig(fname_events, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved plots:\n   {fname_events}")


def main():
    parser = argparse.ArgumentParser(description="Plot diagnostic signal injections")
    parser.add_argument("--output", type=str, required=True, help="Output directory for plots")
    parser.add_argument("--pulseshape", type=str, choices=["mvm", "sine"], required=True)
    parser.add_argument("--method", type=str, choices=["classic", "base"], required=True)
    parser.add_argument("--frequency", type=float, required=True)
    parser.add_argument("--phi", type=float, default=0.0)
    parser.add_argument("--a", type=float, required=True, help="Amplitude")
    parser.add_argument("--baseline", type=float, default=0.0)
    parser.add_argument("--kappa", type=float, default=None, help="MVMD shape parameter")
    parser.add_argument("--bin_time", type=float, default=1.0)
    parser.add_argument("--num_periods", type=int, default=5)
    parser.add_argument("--high_res_dt", type=float, default=0.01)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    plot_signal_examples(
        output_dir=args.output,
        pulseshape=args.pulseshape,
        method=args.method,
        frequency=args.frequency,
        phi=args.phi,
        a=args.a,
        baseline=args.baseline,
        kappa=args.kappa,
        bin_time=args.bin_time,
        num_periods=args.num_periods,
        high_res_dt=args.high_res_dt
    )


if __name__ == "__main__":
    main()
