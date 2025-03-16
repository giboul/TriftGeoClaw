#!/usr/bin/env python
"""
Module for writing the qinit.xyz describing
the initial surface level of the lake in TSUL.

Running
-------
Either call `main(plot: bool=False)`
or run this script with `python makeqinit_tsul.py [--plot]`
"""
from pathlib import Path
from yaml import safe_load
import utils
import numpy as np
from matplotlib import pyplot as plt


projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    TOPM = config["TOPM"]
    TSUL = config["TSUL"]


def main(plot=False):
    # Read the bathymetry
    x, y, Z = utils.read_asc(projdir/TOPM["bathymetry"])
    extent = x.min(), x.max(), y.min(), y.max()
    if plot:
        plt.imshow(Z, extent=extent)
    X, Y = np.meshgrid(x, y, copy=False)

    # Create the lake mask from the seed of fidget around with a GUI
    print("Filling lake", flush=True)
    if 'flood_seed' in TOPM:  # Directly fill 
        seed = (np.abs(y-TOPM['flood_seed'][1]).argmin(),
                np.abs(x-TOPM['flood_seed'][0]).argmin())
        r = int(TOPM["dilation_radius"]/TOPM["resolution"])
        flooded = utils.flood_mask(Z, seed, TOPM["lake_alt"])
        dilated = utils.isotropic_dilation(flooded, r)
    else:  # Interactive fill with mpl figure
        flooded, dilated = utils.pick_seed(Z, x, y, TOPM['resolution'], TOPM['lake_alt'])

    # Save the qinit.xyz
    print(f"Saving {projdir / TSUL['qinit']}", flush=True)
    Z[flooded] = TOPM["lake_alt"]
    Z[~flooded] = 0.
    np.savetxt(projdir / TSUL['qinit'], np.column_stack((
        X.flatten(), Y.flatten(), np.where(dilated, TOPM["lake_alt"], Z.min()).flatten()
    )), fmt="%.9e")

    # Find the contour of the lake for later processing (energy and momentum analysis)
    contour = utils.find_contour(dilated.T, extent, x.size, y.size)
    np.savetxt(projdir / "tsul" / "contour.xy", contour)

    # And write the computational bounds around the lake for fixed grid output
    lake_extent = utils.expand_bounds(
        X[flooded].min(), X[flooded].max(),
        Y[flooded].min(), Y[flooded].max(),
        rel_margin=1/20,
        abs_margin=10
    )
    np.savetxt(projdir / "topm" / "lake_extent.txt", lake_extent)

    if plot:
        plt.imshow(np.where(Z>0, Z, np.nan), extent=extent)
        plt.show()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser() 
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    main(**args.__dict__)

