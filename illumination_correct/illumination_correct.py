'''
Run function for illumination correction of LLAMAS data cubes. 
This script reads a list of input data cubes from a .txt file, 
applies a median filter for illumination correction, and saves 
the corrected data cubes along with the median filtered cubes 
to the specified output directory.
'''

import numpy as np
import matplotlib.pyplot as plt
import argparse
import json
from pathlib import Path
from astropy.io import fits
import warnings
warnings.filterwarnings("ignore")
from copy import deepcopy
from median_filter import median_2dfilter
from treat_cube import treat_cube
from rotation_solve import make_footprint
from rotation_gui import launch_gui

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply illumination correction to LLAMAS data cubes.")
    parser.add_argument("input_files", type=str, help="Path to the .txt list file specifying the input data cubes.")
    parser.add_argument("--params", type=str, help="Path to the .json file with the footprint shape and angle saved by rotation_gui.py. If not given, the GUI is launched first to define the footprint and saves footprint.json next to the input list.")
    parser.add_argument("--pixscale", type=float, default=0.5, help="Spaxel size in arcsec of the cubes built from the RSS files.")
    parser.add_argument("--skip_cube_creation", action="store_true", help="Skip cube creation if the input_files.txt contains paths to already created cubes.")
    parser.add_argument("--output_dir", type=str, help="Path to directory for saving the output files.")
    args = parser.parse_args()

    # Load
    with open(args.input_files, "r") as f:
        input_files = [line.strip() for line in f.readlines()]

    # Footprint of the median filter, as accepted in rotation_gui.py.
    # Without --params, launch the GUI on the first cube to define it
    if args.params:
        params_path = Path(args.params)
    else:
        params_path = Path(args.input_files).with_name("footprint.json")
        first_path = Path(input_files[0])
        gui_output_dir = Path(args.output_dir) if args.output_dir else first_path.parent
        old_mtime = params_path.stat().st_mtime if params_path.exists() else None
        launch_gui(first_path, gui_output_dir, params_path, args.pixscale, args.skip_cube_creation)
        # An old footprint.json must not be reused if the GUI was closed without clicking Accept
        if params_path.exists() and params_path.stat().st_mtime == old_mtime:
            raise SystemExit("The GUI was closed without clicking Accept. Nothing was done.")
    if not params_path.exists():
        raise SystemExit(f"{params_path} not found. Click Accept in the GUI, or pass --params with a file saved by rotation_gui.py.")
    with open(params_path, "r") as f:
        params = json.load(f)
    footprint = make_footprint(tuple(params["shape"]), params["angle"])
    print(f"Footprint shape {tuple(params['shape'])}, angle {params['angle']:.1f} deg, from {params_path}")

    # Process each input file
    for input_file in input_files:
        input_path = Path(input_file)
        if not input_path.exists():
            print(f"Input file {input_file} not found. Skipping.")
            continue

        # Build the data cube from the RSS file if not skipping cube creation
        if not args.skip_cube_creation:
            output_dir = Path(args.output_dir) if args.output_dir else input_path.parent
            cube_path = treat_cube(input_path, output_dir, args.pixscale)
        else:
            cube_path = input_path
            output_dir = Path(args.output_dir) if args.output_dir else input_path.parent

        # Load data cube
        with fits.open(cube_path) as hdul:
            data_cube = hdul[0].data
            header = hdul[0].header
            header1 = deepcopy(header)
            header2 = deepcopy(header)

        print(f"Processing {input_file} with shape {data_cube.shape}")

        # Apply median filter for illumination correction
        median_cube = median_2dfilter(data_cube, footprint)
        # Add header card to indicate median filter was applied
        header1['ILLUM_CORR'] = 'Median filter applied for illumination correction'

        # Save median filtered cube
        # median cube must be subtracted from the original data cube to correct for illumination
        corrected_cube = data_cube - median_cube# Does this work? I think it should be data_cube - median_cube, but let's check the dimensions
        header2['ILLUM_CORR'] = 'Median filtered cube subtracted from original data cube'
        header2['ORIG_CUBE'] = f'Original data cube {input_path.name}'
        median_output_file = output_dir / f"{input_path.stem}_median.fits"
        fits.writeto(median_output_file, median_cube, header1, overwrite=True)

        # Save corrected data cube
        corrected_output_file = output_dir / f"{input_path.stem}_corrected.fits"
        fits.writeto(corrected_output_file, corrected_cube, header2, overwrite=True)
        print(f"Saved corrected data cube to {corrected_output_file}")