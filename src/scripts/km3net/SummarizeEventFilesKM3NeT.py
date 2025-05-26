""" Summarizing CombineEventListsKM3NeT.py output files: *_total_events.txt.

Usage: CombineEventListsKM3NeT.py

Options:
  -h --help                              Show this help message
"""


from docopt import docopt
import os
import numpy as np

import datetime


def summarize_event_files(folder_path, file_suffix="total_events.txt", output_file="summary.txt"):
    summary_data = []

    for filename in os.listdir(folder_path):
        if filename.endswith(file_suffix):
            file_path = os.path.join(folder_path, filename)
            try:
                with open(file_path, 'r') as f:
                    lines = [line.strip() for line in f.readlines()]
                    if len(lines) < 4:
                        print(f"Skipping incomplete file: {filename}")
                        continue

                    n_total = int(lines[0])
                    n_estimated = float(lines[1])
                    dt_str = lines[2]
                    source_name = lines[3]
                    opening_angle = lines[4]
                    detector = lines[5]
                    filestype = lines[6]
                    

                    dt = datetime.datetime.strptime(dt_str, "%Y-%m-%d %H:%M:%S")
                    summary_data.append((dt, n_total, n_estimated, source_name, opening_angle, detector, filestype, filename))

            except Exception as e:
                print(f"Error processing {filename}: {e}")

    # Sort by datetime (descending: newest first)
    summary_data.sort(reverse=True, key=lambda x: x[0])

    # Write summary file
    #with open(os.path.join(folder_path, output_file), 'w') as out:
    #    out.write("# Datetime\t\t\tTotal_Events\tEstimated_Neutrinos\tSource\t\tFilename\n")
    #    for dt, total, estimated, source, fname in summary_data:
    #        out.write(f"{dt}\t\t{total}\t\t{estimated:.2f}\t\t\t{source}\t{fname}\n")

    with open(os.path.join(folder_path, output_file), 'w') as out:
        out.write(f"# {'Datetime':<20} {'Total_Events':>12} {'Estimated_Neutrinos':>20}  {'Source':<15} {'Min_Opening_Angle':<20} {'Detector':<20} {'File_Type':<15} Filename\n")
        for dt, total, estimated, source, opening_angle, detector, filestype, fname in summary_data:
            out.write(f"{dt.strftime('%Y-%m-%d %H:%M:%S')}  {total:>12}  {estimated:>20.2f}  {source:<15} {opening_angle:<20} {detector:<20} {filestype:<15} {fname}\n")

    print(f"Summary written to {output_file}")

def main():
    #arguments = docopt(__doc__)
    #folder_path_read  = arguments['--folder_path']
    folder_path = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists"
    summarize_event_files(folder_path)

if __name__ == "__main__":
    main()