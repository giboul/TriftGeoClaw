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
import topo_utils
import numpy as np
from matplotlib import pyplot as plt


with open("config.yaml") as file:
    config = safe_load(file)


def main(plot=False):
    # Read the bathymetry
    x, y, Z = topo_utils.read_asc(Path(config["bathymetry"]).expanduser())
    resolution = abs(x[1]-x[0])
    extent = x.min(), x.max(), y.min(), y.max()
    if plot:
        plt.imshow(Z, extent=extent)
    X, Y = np.meshgrid(x, y, copy=False)

    # Create the lake mask from the seed or fidget around with a GUI
    print("Filling lake", flush=True)
    if 'flood_seed' in config:  # Directly fill 
        seed = (np.abs(y-config['flood_seed'][1]).argmin(),
                np.abs(x-config['flood_seed'][0]).argmin())
        r = int(config["dilation_radius"]/resolution)
        flooded = topo_utils.flood_mask(Z, seed, config["lake_alt"])
        dilated = topo_utils.isotropic_dilation(flooded, r)
    else:  # Interactive fill with mpl figure
        flooded, dilated = topo_utils.pick_seed(Z, x, y, config['lake_alt'])

    # Save the qinit.xyz
    print(f"Saving qinit.xyz", flush=True)
    Z[flooded] = config["lake_alt"]
    Z[~flooded] = 0.
    np.savetxt('qinit.xyz', np.column_stack((
        X.flatten(), Y.flatten(), np.where(dilated, config["lake_alt"], Z.min()).flatten()
    )), fmt="%.9e")

    # Find the contour of the lake for later processing (energy and momentum analysis)
    contour = topo_utils.find_contour(dilated.T, extent, x.size, y.size)
    np.savetxt("contour.xy", contour)

    if plot:
        plt.imshow(np.where(Z>0, Z, np.nan), extent=extent)
        plt.show()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser() 
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    main(**args.__dict__)

