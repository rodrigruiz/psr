"""
Usage:
  TestScript2.py -i <input_files>... -o <output_dir>

Options:
  -i <input_files>    Input file paths.
  -o <output_dir>     Output directory.
"""

from docopt import docopt  # Ensure we import the function directly

def main():
    import sys
    # Print sys.argv for debugging
    print("sys.argv:", sys.argv)

    # Print docopt version
    import docopt
    print("docopt version:", docopt.__version__)

    # Getting Key-Argument-Pairs that are passed to the script
    arguments = docopt(__doc__, argv=sys.argv[1:])

    # Debug: Print the parsed arguments
    print("Arguments:", arguments)

    input_files = arguments['<input_files>']
    output_dir = arguments['<output_dir>']

    # Debug: Check input files and output directory
    print("Input files:", input_files)
    print("Output directory:", output_dir)

if __name__ == "__main__":
    main()
