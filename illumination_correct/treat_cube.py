'''
Run function to build LLAMAS data cubes from their RSS files.
This script reads a list of input RSS files from a .txt list file, builds
a data cube for each one with the llamas-pyjamas pipeline functions
(RSS -> super-RSS -> cube), and saves the cubes to the specified output
directory. These cubes are the input for the illumination correction.
Requires the llamas-pyjamas environment.
'''

import numpy as np
import argparse
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
from llamas_pyjamas.Combine.superRSS import build_super_rss
from llamas_pyjamas.Combine.cube import combine_cube

def treat_cube(input_path, output_dir, pixscale=0.5):
    """
    Build the data cube of one RSS file and save it to output_dir.
    Returns the path of the saved cube.
    """
    input_path = Path(input_path)

    # Exposure prefix and channel from the file name, e.g. LLAMAS_<date>_RSS_blue.fits
    prefix, channel = input_path.stem.split("_RSS_")

    # Build the cube for this exposure and channel
    super_rss = build_super_rss([str(input_path)], channels=[channel])
    cube = combine_cube(super_rss, channel, pixscale=pixscale)

    # Save cube (SCI primary + VAR, COVERAGE, NEXP, WAVELENGTH extensions)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cube_output_file = output_dir / f"{prefix}_cube_{channel}.fits"
    cube.write(cube_output_file)
    print(f"Saved cube {cube.data.shape} to {cube_output_file}")

    return cube_output_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build LLAMAS data cubes from RSS files.")
    parser.add_argument("input_files", type=str, help="Path to the .txt list file specifying the input RSS files.")
    parser.add_argument("--pixscale", type=float, default=0.5, help="Spaxel size in arcsec.")
    parser.add_argument("--output_dir", type=str, help="Path to directory for saving the output cubes.")
    args = parser.parse_args()

    # Load
    with open(args.input_files, "r") as f:
        input_files = [line.strip() for line in f.readlines()]

    # Process each input file
    for input_file in input_files:
        input_path = Path(input_file)
        if not input_path.exists():
            print(f"Input file {input_file} not found. Skipping.")
            continue

        output_dir = Path(args.output_dir) if args.output_dir else input_path.parent / "cubes"
        treat_cube(input_path, output_dir, args.pixscale)
