'''
GUI to define the shape and rotation of the median filter footprint.
This script displays the median white light image of ONE cube (the first in the
input list), for a wavelength range selected by the user. The footprint box
(shape and rotation angle set by the user) is drawn on top of the image, so it
can be aligned by eye with the slice-looking illumination pattern.
Requires the llamas-pyjamas environment.
'''

import numpy as np
import argparse
import json
from pathlib import Path
import tkinter as tk
from tkinter import ttk
from astropy.io import fits
from matplotlib.figure import Figure
from matplotlib.colors import ListedColormap
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import warnings
warnings.filterwarnings("ignore")
import rotation_solve
from rotation_solve import whitelight, make_footprint

def load_cube(cube_path):
    """
    Load a cube saved by treat_cube.py.
    Returns the (nwave, ny, nx) data and the wavelength axis in Angstrom.
    """
    with fits.open(cube_path) as hdul:
        data_cube = hdul[0].data
        # CRVAL3/CDELT3 are rescaled by astropy, the WAVELENGTH table is the reliable axis
        wave = np.asarray(hdul["WAVELENGTH"].data["WAVELENGTH"], dtype=float)
    return data_cube, wave

class RotationGUI:
    def __init__(self, root, data_cube, wave, title="", params_path="footprint.json"):
        self.root = root
        self.params_path = Path(params_path)
        self.cube_name = title
        self.data_cube = data_cube
        self.wave = wave
        root.title(f"Footprint rotation {title}")

        # Variables set by the widgets
        self.wmin = tk.StringVar(value=f"{wave[0]:.1f}")
        self.wmax = tk.StringVar(value=f"{wave[-1]:.1f}")
        self.nx = tk.StringVar(value="3")
        self.ny = tk.StringVar(value="101")
        self.angle = tk.DoubleVar(value=0.0)
        self.status = tk.StringVar(value="Click on the image to move the box.")

        self.build_widgets()

        # Footprint center in image pixels (x, y), default: middle of the image
        ny_img, nx_img = data_cube.shape[1:]
        self.center = ((nx_img - 1) / 2, (ny_img - 1) / 2)

        self.update_image()

    def build_widgets(self):
        controls = ttk.Frame(self.root, padding=8)
        controls.grid(row=0, column=0, sticky="n")

        # Wavelength range for the white light image
        ttk.Label(controls, text="Wavelength range (A)").grid(row=0, column=0, columnspan=2, sticky="w")
        ttk.Entry(controls, textvariable=self.wmin, width=9).grid(row=1, column=0)
        ttk.Entry(controls, textvariable=self.wmax, width=9).grid(row=1, column=1)
        ttk.Button(controls, text="Update image", command=self.update_image).grid(row=2, column=0, columnspan=2, pady=4)

        # Footprint shape (x, y)
        ttk.Label(controls, text="Footprint shape (x, y)").grid(row=3, column=0, columnspan=2, sticky="w", pady=(10, 0))
        entry_x = ttk.Entry(controls, textvariable=self.nx, width=9)
        entry_y = ttk.Entry(controls, textvariable=self.ny, width=9)
        entry_x.grid(row=4, column=0)
        entry_y.grid(row=4, column=1)
        entry_x.bind("<Return>", self.draw)
        entry_y.bind("<Return>", self.draw)

        # Rotation angle: slider and entry share the same variable
        ttk.Label(controls, text="Rotation angle (deg, counterclockwise)").grid(row=5, column=0, columnspan=2, sticky="w", pady=(10, 0))
        tk.Scale(controls, from_=0, to=180, resolution=0.1, orient="horizontal", length=180,
                 variable=self.angle, command=self.draw).grid(row=6, column=0, columnspan=2)
        entry_angle = ttk.Entry(controls, textvariable=self.angle, width=9)
        entry_angle.grid(row=7, column=0, columnspan=2)
        entry_angle.bind("<Return>", self.draw)

        ttk.Button(controls, text="Solve Rotation", command=self.solve).grid(row=8, column=0, columnspan=2, pady=(12, 4))
        ttk.Button(controls, text="Accept", command=self.accept).grid(row=9, column=0, columnspan=2, pady=4)
        ttk.Label(controls, textvariable=self.status, wraplength=190).grid(row=10, column=0, columnspan=2, sticky="w", pady=(8, 0))

        # Image
        self.fig = Figure(figsize=(6, 6))
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().grid(row=0, column=1)
        self.canvas.mpl_connect("button_press_event", self.on_click)

    def update_image(self, *args):
        """Recompute the white light image for the selected wavelength range."""
        try:
            wmin, wmax = float(self.wmin.get()), float(self.wmax.get())
        except ValueError:
            self.status.set("Wavelength range must be numbers.")
            return
        # Angstrom -> plane indices
        i0 = np.searchsorted(self.wave, wmin)
        i1 = np.searchsorted(self.wave, wmax, side="right")
        if i1 - i0 < 1:
            self.status.set("Empty wavelength range.")
            return
        self.image = whitelight(self.data_cube, i0, i1)
        self.status.set(f"White light of {i1 - i0} planes ({self.wave[i0]:.0f}-{self.wave[i1 - 1]:.0f} A).")
        self.draw()

    def draw(self, *args):
        """Redraw the white light image with the footprint on top."""
        try:
            shape = (int(self.nx.get()), int(self.ny.get()))
            angle = self.angle.get()
        except (ValueError, tk.TclError):
            self.status.set("Shape must be integers and angle a number.")
            return
        if min(shape) < 1:
            self.status.set("Shape must be at least 1.")
            return

        footprint = make_footprint(shape, angle)
        side = footprint.shape[0]
        cx, cy = self.center
        ny_img, nx_img = self.image.shape

        self.ax.clear()
        vmin, vmax = np.nanpercentile(self.image, [2, 98])
        self.ax.imshow(self.image, origin="lower", vmin=vmin, vmax=vmax)
        # Footprint pixels in red, centered on (cx, cy); pixel i covers [i - 0.5, i + 0.5]
        self.ax.imshow(np.ma.masked_where(~footprint, footprint), origin="lower",
                       cmap=ListedColormap(["red"]), alpha=0.7,
                       extent=(cx - side / 2, cx + side / 2, cy - side / 2, cy + side / 2))
        self.ax.set_xlim(-0.5, nx_img - 0.5)
        self.ax.set_ylim(-0.5, ny_img - 0.5)
        self.ax.set_title(f"shape (x, y) = {shape}, angle = {angle:.1f} deg")
        self.canvas.draw()

    def on_click(self, event):
        """Move the footprint to the clicked position."""
        if event.inaxes is self.ax:
            self.center = (event.xdata, event.ydata)
            self.draw()

    def solve(self):
        """Fit the rotation angle to the white light image shown and set the slider to it."""
        try:
            shape = (int(self.nx.get()), int(self.ny.get()))
        except ValueError:
            self.status.set("Shape must be integers.")
            return
        if min(shape) < 1:
            self.status.set("Shape must be at least 1.")
            return

        # The fit takes several seconds and blocks the window, so show a message first
        self.status.set("Solving rotation...")
        self.root.update_idletasks()
        angle, contrast = rotation_solve.solve_rotation(self.image, shape)

        self.angle.set(round(angle, 1))
        self.draw()
        message = f"Solved angle = {angle:.1f} deg (contrast {contrast:.2f})."
        if contrast < 0.05:
            message += " Low contrast: no clear pattern in this wavelength range."
        self.status.set(message)

    def accept(self):
        """Save the shape and angle to the params JSON file read by illumination_correct.py, and close."""
        try:
            shape = [int(self.nx.get()), int(self.ny.get())]
            angle = self.angle.get()
        except (ValueError, tk.TclError):
            self.status.set("Shape must be integers and angle a number.")
            return
        if min(shape) < 1:
            self.status.set("Shape must be at least 1.")
            return

        params = {"shape": shape, "angle": angle,
                  "wmin": float(self.wmin.get()), "wmax": float(self.wmax.get()),
                  "cube": self.cube_name}
        with open(self.params_path, "w") as f:
            json.dump(params, f, indent=2)
        print(f"Saved shape {tuple(shape)} and angle {angle:.1f} deg to {self.params_path}")
        self.root.destroy()

def launch_gui(input_file, output_dir, params_path, pixscale=0.5, skip_cube_creation=False):
    """
    Open the GUI on the cube of one input file and wait until the window is closed.
    The cube is built from the RSS file (reusing it if already in output_dir), or the
    input file is used directly as the cube if skip_cube_creation is set.
    The shape and angle accepted by the user are saved to params_path.
    """
    input_path = Path(input_file)
    output_dir = Path(output_dir)
    if skip_cube_creation:
        cube_path = input_path
    else:
        prefix, channel = input_path.stem.split("_RSS_")
        cube_path = output_dir / f"{prefix}_cube_{channel}.fits"
        if not cube_path.exists():
            from treat_cube import treat_cube
            cube_path = treat_cube(input_path, output_dir, pixscale)

    data_cube, wave = load_cube(cube_path)
    print(f"Loaded {cube_path.name} with shape {data_cube.shape}")

    root = tk.Tk()
    RotationGUI(root, data_cube, wave, title=cube_path.name, params_path=params_path)
    root.mainloop()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GUI to set the shape and rotation of the median filter footprint.")
    parser.add_argument("input_files", type=str, help="Path to the .txt list file specifying the input RSS files. The first one is used.")
    parser.add_argument("--pixscale", type=float, default=0.5, help="Spaxel size in arcsec of the cube built from the RSS file.")
    parser.add_argument("--skip_cube_creation", action="store_true", help="Skip cube creation if the input_files.txt contains paths to already created cubes.")
    parser.add_argument("--output_dir", type=str, help="Path to directory where the cubes are saved (as in illumination_correct.py).")
    parser.add_argument("--params", type=str, help="Path of the .json file to save the footprint shape and angle to (default: footprint.json next to the input list).")
    args = parser.parse_args()

    # Load
    with open(args.input_files, "r") as f:
        input_files = [line.strip() for line in f.readlines()]

    # Use the first file
    input_path = Path(input_files[0])
    output_dir = Path(args.output_dir) if args.output_dir else input_path.parent
    params_path = Path(args.params) if args.params else Path(args.input_files).with_name("footprint.json")
    launch_gui(input_path, output_dir, params_path, args.pixscale, args.skip_cube_creation)
