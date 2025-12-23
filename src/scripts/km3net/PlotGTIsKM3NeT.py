"""
Plot lightcurve with GTIs overlay.

Usage:
    PlotGTIsKM3NeT.py -e EVENTLIST -g GTIFILE -o OUTPUT_PNG [--sourcename=<sourcename>] [--df=<df>] 

Options:
    -h --help                Show this screen
    -e --eventlist EVENTLIST Path to combined eventlist HDF5 file
    -g --gtifile GTIFILE     Path to GTIs file created by saveGTIs
    -o --output OUTPUT_PNG   Output PNG filename
       --sourcename=<sourcename>     Source name [default: Vela]
       --df=<float>                  Time bin size in seconds [default: 60]
"""

from docopt import docopt
import h5py
import numpy as np
import matplotlib.pyplot as plt
from stingray import Lightcurve
import plens.EventList as EL
from epochfolding.gtis import loadGTIs

def main():
    args = docopt(__doc__)
    eventlist_file = str(args["--eventlist"])
    gtifile = str(args["--gtifile"])
    output_png = str(args["--output"])
    sourcename = str(args["--sourcename"])
    df = float(args["--df"])

    # Load event list
    with h5py.File(eventlist_file, "r") as f:
        eventlist = EL.readEventList(f)

    if len(eventlist) == 0:
        print("No events found in the eventlist file.")
        return

    # Load GTIs (assumes saved with numpy.savetxt or similar)
    gtis = loadGTIs(gtifile)
    #if gtis.ndim == 1 and gtis.size == 2:  # Single GTI
    #    gtis = gtis[np.newaxis, :]

    # Make lightcurve
    lc = Lightcurve.make_lightcurve(eventlist["time"].value, dt=df)

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(lc.time, lc.counts, label="Counts", color="darkblue")

    duration_list = []

    for start, end in gtis:
        plt.axvline(start, linestyle="dotted", color="gray", alpha=0.8)
        plt.axvline(end, linestyle="dotted", color="gray", alpha=0.8)
        plt.axvspan(start, end, color="gray", alpha=0.2)
        duration = end - start
        print("Duration: ", duration)
        duration_list.append(duration)

    active_time = np.sum(np.array(duration_list))
    active_time_d = active_time/(3600*24)
    print(f"active_time : {active_time:.2f}s - {active_time_d:.2f} d")
    plt.xlabel("Time (s)", fontsize=14)
    plt.ylabel("Counts", fontsize=14)
    plt.title(f"Lightcurve with highlighted GTIs (ARCA: {sourcename}, 30deg)", fontsize=16)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_png)
    plt.close()

    print(f"Saved GTI plot to {output_png}")


if __name__ == "__main__":
    main()
