""" Plots Sky Maps of the Input EventLists an highlights the source from the Source Specs File.

Usage: PlotSkyMapsKM3NeT.py -i INPUT_FILE -o OUTPUT_DIR -s SOURCE_SPECS_FILE [--radius=<radius>] [--num_bins=<num_bins>] [--detectorname=<detectorname>]

Options:
  -h --help                              Show this help message
  -i --input_file INPUT_FILE             Input file
  -o --output_dir OUTPUT_DIR             Output directory  
  -s --source SOURCE_SPECS_FILE          Hdf5 file containing information about the source of interest (ra, dec, P_orb, ...) 
     --radius=<float>                    Search Cone Minimal Radius [default: 5.]
     --num_bins=<int>                    Number of Bins along the SkyMap axes for the histogram [default: 200]
     --detectorname=<detectorname>       Detectorname ('arca' or 'orca) [default: arca]
  """ 

# python3 psr/src/scripts/CombineEventListsKM3NeT.py -i "/home/hpc/capn/capn107h/software/correctedeventlistTestOutput/*" -o combinedeventlistTestOutput/


from docopt import docopt
import os
import glob
import h5py
import numpy as np
from astropy.table import vstack, Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.units import degree
import plens.EventList as EL
import pylab as pl

import km3astro.plot as km3plt
from km3astro.coord import local_event
import matplotlib.colors as mcolors
from scipy.ndimage import gaussian_filter

from stingray import EventList, Lightcurve
import warnings
import matplotlib.pyplot as plt
from astropy.time import Time
from h5py import File

def plot_equatorial_custom(
    evts,
    projection="aitoff",
    ax=None,
    marker="o",
    alpha=0.8,
    markersize =2,
    color='darkblue',
    no_face_colors = False,
    adjust_subplots=True,
    **kwargs,
):
    if no_face_colors is True: 
        kwargs['facecolors']='none'
        kwargs['edgecolors']=color
        kwargs['s']=markersize*12
    else: 
        kwargs['c']=color
        kwargs['s']=markersize

    ra, dec = km3plt.ra_dec(evts)
    if ax is None:
        _, ax = km3plt.projection_axes(projection=projection)
    ax.scatter(np.array(ra), np.array(dec), marker=marker, alpha=alpha, **kwargs)
    if adjust_subplots:
        plt.subplots_adjust(top=0.95, bottom=0.0)
    return ax

def plot_equatorial_skymap(theta,phi,event_times,energy,source_ra,source_dec,source_name, output_file, title_prefix="", radius=5*u.deg, circle=False,energy_max_size = None, detectorname = 'arca'):

    source_coord = SkyCoord(ra=source_ra * u.deg, dec=source_dec * u.deg)

    
    total_events = len(event_times)
    
    # Convert detector frame angles to sky coordinates (ICRS)
    #event_location_icrs = local_event(theta, phi, astropy_event_times, 'arca').transform_to('icrs')
    event_location_icrs = local_event(theta.value, phi.value, event_times, detectorname).transform_to('icrs')

    # Sort events by energy (low to high)
    sorted_indices = np.argsort(energy)  # Sort energy in ascending order
    event_location_icrs = event_location_icrs[sorted_indices]
    energy = energy[sorted_indices]

    distances = event_location_icrs.separation(source_coord).degree
    
    # Count selected events
    count_in_circle = np.sum(distances <= radius.value)
    selection_ratio = count_in_circle / total_events if total_events > 0 else 0

    # Normalize energy for color coding
    energy_min = energy.min()
    energy_max = energy.max()
    #if energy_max is None: energy_max = energy.max()
    energy_normalized = (energy - energy_min) / (energy_max - energy_min)
    # Ensure values above max_energy are set to 1
    #energy_normalized = np.where(energy > energy_max, 1, energy_normalized)

    print(energy_max)

    # Map energy values to sizes
    size_min, size_max = 5, 30
    energy_size_normalized = energy_normalized #(energy - energy_min) / (energy_max_size - energy_min)
    sizes = size_min + energy_size_normalized * (size_max - size_min)
    sizes = np.where(sizes > size_max,size_max, sizes)

    # Create RGBA colors
    cmap = plt.cm.viridis
    colors = cmap(energy_normalized)
    #colors[:, -1] = 0.2 + 0.8 * (energy_normalized)  # Alpha increases with energy

    #sizes = sizes[sorted_indices]


    fig, ax = km3plt.projection_axes(projection="aitoff", figsize=(10, 5))
    ax = plot_equatorial_custom(event_location_icrs,ax=ax,markersize=sizes,edgecolors="none",color=colors,label="Events")
    plot_equatorial_custom(source_coord, markersize=5, color='deeppink',marker="o", ax=ax, label=source_name, no_face_colors=True)
    #ax.legend()
    

    # Colorbar setup
    norm_log = mcolors.LogNorm(vmin=max(energy_min, 1e-5), vmax=energy_max)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm_log)
    cbar = plt.colorbar(sm, ax=ax, orientation="horizontal", pad=0.07)
    cbar.set_label("Neutrino Energy")

    if circle is True:
        circle_radius = np.radians(radius.to(u.deg).value)
        ax.add_patch(plt.Circle(
            (np.radians(source_coord.ra.wrap_at(180 * u.degree).degree), np.radians(source_coord.dec.degree)),
            circle_radius,
            transform=ax.transData,
            edgecolor="deeppink",
            facecolor="none",
            linestyle="--",
            linewidth=2,
            alpha=1,
            zorder=2,
        ))
    
        # Create an additional axis for legend and text box
        right_ax = fig.add_axes([0.82, 0.3, 0.15, 0.4])  # x, y, width, height
        right_ax.axis("off")  # Hide the axis



        # Add the text box right below the legend
        text_box_content = (f"Total events: {total_events}\n"
                            f"Events within {radius.value:.1f}°: {count_in_circle}\n"
                            f"Selection ratio: {selection_ratio*100:.2f} %")
        right_ax.text(
            0.5, 0.7, text_box_content, fontsize=10, transform=right_ax.transAxes,
            verticalalignment="top", horizontalalignment="center",
            bbox=dict(facecolor="white", edgecolor="black", alpha=0.9)
        )
            # Add the legend at the top of the right panel

        legend_handles, legend_labels = ax.get_legend_handles_labels()
        legend = right_ax.legend(
            legend_handles, legend_labels, loc="upper center", frameon=True, fontsize=10
        )
    else: ax.legend(loc='upper right')

    ax.set_title(f"{title_prefix} Sky Map for {source_name} \n {total_events} Total Events - {radius} Search Cone", fontsize=16, pad=30)
    
    output_plot = output_file + "_SkyMap.png"
    plt.savefig(output_plot, bbox_inches='tight')
    print(f"Plot saved: {output_plot}")


def plot_equatorial_hist_skymap(theta,phi,event_times,num_bins,source_ra,source_dec,source_name, output_file, title_prefix="", radius=5*u.deg, circle=False,energy_max_size = None, zero_white_bins = True, smooth_hist = True, detectorname='arca'):

    source_coord = SkyCoord(ra=source_ra * u.deg, dec=source_dec * u.deg)

    
    total_events = len(event_times)
    
    # Convert detector frame angles to sky coordinates (ICRS)
    #event_location_icrs = local_event(theta, phi, astropy_event_times, 'arca').transform_to('icrs')
    event_location_icrs = local_event(theta.value, phi.value, event_times, detectorname).transform_to('icrs')

    distances = event_location_icrs.separation(source_coord).degree
    
    # Count selected events
    count_in_circle = np.sum(distances <= radius.value)
    selection_ratio = count_in_circle / total_events if total_events > 0 else 0
    ra_vals, dec_vals = km3plt.ra_dec(event_location_icrs)

    # Define bin edges
    ra_bins = np.linspace(-np.pi, np.pi, num_bins + 1)  # RA in [-π, π]
    dec_bins = np.linspace(-np.pi/2, np.pi/2, num_bins + 1)  # Dec in [-π/2, π/2]

    # Compute 2D histogram
    hist, ra_edges, dec_edges = np.histogram2d(ra_vals, dec_vals, bins=[ra_bins, dec_bins])

    if smooth_hist is True:
        hist_g = gaussian_filter(hist, sigma=1.0)
        if zero_white_bins is True: hist_plot = np.ma.masked_where(hist == 0, hist_g)
        else: hist_plot = hist_g
    else:
        if zero_white_bins is True: hist_plot = np.ma.masked_where(hist == 0, hist)
        else: hist_plot = hist

    # Convert bin centers
    ra_centers = 0.5 * (ra_edges[:-1] + ra_edges[1:])
    dec_centers = 0.5 * (dec_edges[:-1] + dec_edges[1:])
    ra_grid, dec_grid = np.meshgrid(ra_centers, dec_centers, indexing="ij")


    fig, ax = km3plt.projection_axes(projection="aitoff", figsize=(10, 5))

    pc = ax.pcolormesh(ra_bins[:-1], dec_bins[:-1], hist_plot.T, shading="auto")
    plot_equatorial_custom(source_coord, markersize=5, color='deeppink',marker="o", ax=ax, label=source_name, no_face_colors=True)
    
    cbar = plt.colorbar(pc, ax=ax, orientation="horizontal", pad=0.07)
    cbar.set_label("Counts per Bin")
    ax.legend(loc='upper right')

    ax.set_title(f"{title_prefix} Histogram Sky Map for {source_name} \n {total_events} Total Events - {radius} Search Cone - {num_bins} Bins", fontsize=16, pad=30)
    
    # maybe change something here?
    if smooth_hist is True:
        output_plot = output_file + "_" + str(num_bins) + "bins" + "_smooth_HistoSkyMap.png"
    else:
        output_plot = output_file + "_" + str(num_bins) + "bins" + "_HistoSkyMap.png"
    plt.savefig(output_plot, bbox_inches='tight')
    print(f"Plot saved: {output_plot}")

def plot_ra_dec_hist(theta,phi,event_times,num_bins,source_ra,source_dec,source_name, output_file, title_prefix="", radius=5*u.deg, circle=False,energy_max_size = None, zero_white_bins = True, smooth_hist = True, detectorname='arca'):

    source_coord = SkyCoord(ra=source_ra * u.deg, dec=source_dec * u.deg)
    total_events = len(event_times)
    
    # Convert detector frame angles to sky coordinates (ICRS)
    #event_location_icrs = local_event(theta, phi, astropy_event_times, 'arca').transform_to('icrs')
    event_location_icrs = local_event(theta.value, phi.value, event_times, detectorname).transform_to('icrs')

    distances = event_location_icrs.separation(source_coord).degree
    
    # Count selected events
    count_in_circle = np.sum(distances <= radius.value)
    selection_ratio = count_in_circle / total_events if total_events > 0 else 0
    ra_vals, dec_vals = km3plt.ra_dec(event_location_icrs)

    plt.figure(figsize=(8, 5))

    # 2D histogram
    h = plt.hist2d(ra_vals, dec_vals, bins=num_bins, cmap='viridis')
    plt.colorbar(h[3], label='Counts per bin')

    plt.scatter(source_coord.ra.rad, source_coord.dec.rad, s=30, facecolors='none', edgecolors='deeppink', marker='o', label=source_name)

    plt.xlabel("Right Ascension / rad")
    plt.ylabel("Declination / rad")
    plt.title(f"{title_prefix} RA/Dec Histogram - {source_name}")
    if smooth_hist is True:
        output_plot = output_file + f"_{num_bins}bins_smooth_RA_Dec_2DHist.png"
    else:
        output_plot = output_file + f"_{num_bins}bins_RA_Dec_2DHist.png"
    plt.savefig(output_plot, bbox_inches='tight')
    plt.close()

    print(f"2D Histogram plot saved: {output_plot}")

def plot_ra_projection(theta, phi, event_times, num_bins, output_file, source_ra, source_dec, source_name, title_prefix="", detectorname='arca'):
    source_coord = SkyCoord(ra=source_ra * u.deg, dec=source_dec * u.deg)
    
    # Convert to sky coordinates
    event_location_icrs = local_event(theta.value, phi.value, event_times, detectorname).transform_to('icrs')
    ra_vals, dec_vals = km3plt.ra_dec(event_location_icrs)

    # 2D histogram
    hist2d, xedges, yedges = np.histogram2d(ra_vals, dec_vals, bins=num_bins)

    # Project onto RA axis
    ra_projection = hist2d.sum(axis=1)
    ra_bin_centers = 0.5 * (xedges[1:] + xedges[:-1])

    # Plot
    plt.figure(figsize=(8, 4))
    plt.bar(ra_bin_centers, ra_projection, width=(xedges[1] - xedges[0]), color='skyblue')
    plt.xlabel("Right Ascension / rad")
    plt.ylabel("Counts (summed over Declination)")
    plt.title(f"{title_prefix} RA Projection - {source_name}")
    
    output_plot = output_file + f"_{num_bins}bins_RAProjection.png"
    plt.savefig(output_plot, bbox_inches='tight')
    plt.close()

    print(f"RA projection plot saved: {output_plot}")


def plot_dec_projection(theta, phi, event_times, num_bins, output_file, source_ra, source_dec, source_name, title_prefix="", detectorname='arca'):
    source_coord = SkyCoord(ra=source_ra * u.deg, dec=source_dec * u.deg)
    
    # Convert to sky coordinates
    event_location_icrs = local_event(theta.value, phi.value, event_times, detectorname).transform_to('icrs')
    ra_vals, dec_vals = km3plt.ra_dec(event_location_icrs)

    # 2D histogram
    hist2d, xedges, yedges = np.histogram2d(ra_vals, dec_vals, bins=num_bins)

    # Project onto Dec axis
    dec_projection = hist2d.sum(axis=0)
    dec_bin_centers = 0.5 * (yedges[1:] + yedges[:-1])

    # Plot
    plt.figure(figsize=(8, 4))
    plt.bar(dec_bin_centers, dec_projection, width=(yedges[1] - yedges[0]), color='darkblue')
    plt.xlabel("Declination / rad")
    plt.ylabel("Counts (summed over Right Ascension)")
    plt.title(f"{title_prefix} Dec Projection - {source_name}")
    
    output_plot = output_file + f"_{num_bins}bins_DecProjection.png"
    plt.savefig(output_plot, bbox_inches='tight')
    plt.close()

    print(f"Dec projection plot saved: {output_plot}")



def plot_theta_phi_hist(theta,phi,event_times,num_bins,source_ra,source_dec,source_name, output_file, title_prefix="", radius=5*u.deg, circle=False,energy_max_size = None, zero_white_bins = True, smooth_hist = True):



    plt.figure(figsize=(8, 5))

    # 2D histogram
    h = plt.hist2d(phi, theta, bins=num_bins, cmap='viridis')
    plt.colorbar(h[3], label='Counts per bin')

    # Mark the source
    #plt.scatter(source_coord.ra.rad, source_coord.dec.rad, s=30, color='deeppink', marker="o", label=source_name)

    plt.xlabel("Phi")
    plt.ylabel("Theta")
    plt.title(f"{title_prefix} Theta/Phi Histogram - {source_name}")
    if smooth_hist is True:
        output_plot = output_file + f"_{num_bins}bins_smooth_Theta_Phi_2DHist.png"
    else:
        output_plot = output_file + f"_{num_bins}bins_Theta_Phi_2DHist.png"
    plt.savefig(output_plot, bbox_inches='tight')
    plt.close()

    print(f"2D Histogram Theta Phi plot saved: {output_plot}")

import numpy as np
import matplotlib.pyplot as plt

def plot_phi_projection(theta, phi, num_bins, output_file, title_prefix="", source_name=""):
    # Compute 2D histogram (same as original)
    hist2d, xedges, yedges = np.histogram2d(phi, theta, bins=num_bins)

    # Project onto phi axis by summing over theta (axis=1 since theta is along the y-axis)
    phi_projection = hist2d.sum(axis=1)

    # Compute phi bin centers
    phi_bin_centers = 0.5 * (xedges[1:] + xedges[:-1])

    # Plot the projection
    plt.figure(figsize=(8, 4))
    plt.bar(phi_bin_centers, phi_projection, width=(xedges[1]-xedges[0]), color='skyblue')
    plt.xlabel("Phi")
    plt.ylabel("Counts (summed over Theta)")
    plt.title(f"{title_prefix} Phi Projection - {source_name}")
    
    output_plot = output_file + f"_{num_bins}bins_PhiProjection.png"
    plt.savefig(output_plot, bbox_inches='tight')
    plt.close()

    print(f"Phi projection plot saved: {output_plot}")

def plot_theta_projection(theta, phi, num_bins, output_file, title_prefix="", source_name=""):
    # Compute 2D histogram
    hist2d, xedges, yedges = np.histogram2d(phi, theta, bins=num_bins)

    # Project onto theta axis by summing over phi (axis=0 since phi is along the x-axis)
    theta_projection = hist2d.sum(axis=0)

    # Compute theta bin centers
    theta_bin_centers = 0.5 * (yedges[1:] + yedges[:-1])

    # Plot the projection
    plt.figure(figsize=(8, 4))
    plt.bar(theta_bin_centers, theta_projection, width=(yedges[1]-yedges[0]), color='darkblue')
    plt.xlabel("Theta")
    plt.ylabel("Counts (summed over Phi)")
    plt.title(f"{title_prefix} Theta Projection - {source_name}")
    
    output_plot = output_file + f"_{num_bins}bins_ThetaProjection.png"
    plt.savefig(output_plot, bbox_inches='tight')
    plt.close()

    print(f"Theta projection plot saved: {output_plot}")



def plot_ra_dec(theta,phi,event_times,num_bins,source_ra,source_dec,source_name, output_file, title_prefix="", radius=5*u.deg, circle=False,energy_max_size = None, zero_white_bins = True, smooth_hist = True,detectorname='arca'):

    source_coord = SkyCoord(ra=source_ra * u.deg, dec=source_dec * u.deg)
    total_events = len(event_times)
    
    # Convert detector frame angles to sky coordinates (ICRS)
    #event_location_icrs = local_event(theta, phi, astropy_event_times, 'arca').transform_to('icrs')
    event_location_icrs = local_event(theta.value, phi.value, event_times, detectorname).transform_to('icrs')

    distances = event_location_icrs.separation(source_coord).degree
    
    # Count selected events
    count_in_circle = np.sum(distances <= radius.value)
    selection_ratio = count_in_circle / total_events if total_events > 0 else 0
    ra_vals, dec_vals = km3plt.ra_dec(event_location_icrs)

    plt.figure(figsize=(8,5))

    plt.scatter(ra_vals,dec_vals)
    plt.scatter(source_coord.ra.rad,source_coord.dec.rad, s=5, color='deeppink',marker="o", label=source_name)
    plt.xlabel("Right Ascension / rad")
    plt.ylabel("Declination / rad")


    # maybe change something here?
    output_plot = output_file + "_" + str(num_bins) + "bins" + "_Rad_Dec_Plot.png"
    plt.savefig(output_plot, bbox_inches='tight')
    print(f"Plot saved: {output_plot}")

def plot_energy_histogram(energy, nbins, source_name, output_file, logx = True, logy = True):

    print("Energy: ", energy)

    if isinstance(energy, u.Quantity):
        valid_energy = energy[energy > 0 * u.GeV]
    else:
        valid_energy = energy[energy > 0]  # or just energy[energy.value > 0] if it's not Quantity
    if len(valid_energy) == 0:
        raise ValueError("No positive energy values for histogram.")
    print("Valid Energy: ", valid_energy)

    min_log = np.log10(valid_energy.value.min())
    max_log = np.log10(valid_energy.value.max())
    logbins = np.logspace(min_log, max_log, nbins)

    print("Min energy:", valid_energy.min(), "Max energy:", valid_energy.max())



    plt.figure(figsize=(8,5))

    if logx and logy:
        bins = logbins
        plt.hist(energy.value,bins=bins,log=True,color="darkblue")
        pl.gca().set_xscale("log")
    elif logx and not logy:
        bins = logbins
        plt.hist(energy.value,bins=bins,log=False,color="darkblue")
        pl.gca().set_xscale("log")
    elif not logx and logy:
        bins = nbins
        plt.hist(energy.value,bins=bins,log=True,color="darkblue")
    elif not logx and not logy: 
        bins = nbins
        plt.hist(energy.value,bins=bins,log=False,color="darkblue")

    print("Bins: ", bins)

    plt.xlabel("Energy (GeV)")
    plt.ylabel("Counts")
    plt.title(f"Energy of Events around {source_name} ({len(energy)} Events)")

    # maybe change something here?
    output_plot = output_file + "_" + str(nbins) + "bins" + "_EnergyHistogram.png"
    plt.savefig(output_plot, bbox_inches='tight')
    print(f"Plot saved: {output_plot}")

def main():
    arguments = docopt(__doc__)

    input_file = arguments['--input_file']
    output_dir = arguments['--output_dir']

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    with File(arguments['--source'], 'r') as f:
        source_name = f['source_name'][()].decode('utf-8')
        source_ra = f['skycoord/ra'][()]
        source_dec = f['skycoord/dec'][()]
        
    with File(input_file, 'r') as h5_file:
        EventList = Table.read(h5_file)
        
        theta = EventList['theta_detectorframe']
        phi = EventList['phi_detectorframe']
        event_times = EventList['time']
        energy = EventList['energy']
        #detectorname = EventList['detector']
        event_count = len(event_times)

    folder_path, file_name = os.path.split(input_file)
    file_name = os.path.splitext(file_name)[0]
    output_file = f"{output_dir}{file_name}"
    num_bins = int(arguments['--num_bins'])
    radius = float(arguments['--radius'])*u.deg
    detectorname = str(arguments['--detectorname'])


    plot_equatorial_skymap(theta,phi,event_times,energy,source_ra,source_dec,source_name,output_file,radius=radius,detectorname=detectorname)
    plot_equatorial_hist_skymap(theta,phi,event_times,num_bins,source_ra,source_dec,source_name,output_file,radius=radius,detectorname=detectorname) #,smooth_hist=False)
    plot_energy_histogram(energy, 200, source_name, output_file, logx = True, logy = True)
    plot_ra_dec_hist(theta,phi,event_times,num_bins,source_ra,source_dec,source_name,output_file,radius=radius, smooth_hist=False,detectorname=detectorname)
    plot_theta_phi_hist(theta,phi,event_times,num_bins,source_ra,source_dec,source_name,output_file,radius=radius, smooth_hist=False)
    plot_ra_dec(theta,phi,event_times,num_bins,source_ra,source_dec,source_name,output_file,radius=radius,detectorname=detectorname)
    plot_phi_projection(theta, phi, num_bins, output_file,source_name=source_name)
    plot_theta_projection(theta, phi, num_bins, output_file,source_name=source_name)
    plot_ra_projection(theta, phi, event_times, num_bins, output_file, source_ra,source_dec, source_name)
    plot_dec_projection(theta, phi, event_times, num_bins, output_file, source_ra,source_dec, source_name)
    
if __name__ == "__main__":
    main()
    