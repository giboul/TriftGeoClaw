#!/usr/bin/env python
"""
Module for writing the qinit.xyz describing
the initial surface level of the lake.
"""
from typing import Tuple
from pathlib import Path
import topo_utils
import numpy as np
from numpy.typing import ArrayLike
from matplotlib import pyplot as plt


def make_qinit(x: ArrayLike,
               y: ArrayLike,
               bathy: ArrayLike,
               lake_alt: float,
               flood_seed: Tuple[int, int]=None,
               dilation_radius: int=0):

    Z = bathy.copy()

    # Read the bathymetry
    resolution = abs(x[1]-x[0])

    # Create the lake mask from the seed or fidget around with a GUI
    if np.isfinite(flood_seed).all():  # Directly fill
        print("Filling lake", flush=True)
        seed = (np.abs(y-flood_seed[1]+resolution/2).argmin(),
                np.abs(x-flood_seed[0]+resolution/2).argmin())
        r = int(dilation_radius/resolution)
        flooded = topo_utils.flood_mask(Z, seed, lake_alt)
        dilated = topo_utils.isotropic_dilation(flooded, r)
    else:  # Interactive fill with mpl figure
        flooded, dilated = topo_utils.pick_seed(Z, x, y, lake_alt)

    # Save the qinit.xyz
    Z[flooded] = lake_alt
    Z[~flooded] = 0.

    qinit_extent = (
        x[dilated.any(axis=0)].min(), x[dilated.any(axis=0)].max(),
        y[dilated.any(axis=1)].min(), y[dilated.any(axis=1)].max()
    )

    return x, y, np.where(dilated, lake_alt, Z.min()), qinit_extent

def save_qinit(x, y, Z):
    X, Y = np.meshgrid(x, y, indexing="ij")
    print("Saving qinit.xyz", end=" ... ", flush=True)
    np.savetxt('qinit.xyz', np.column_stack((
        X.flatten(), Y.flatten(), Z.flatten()
    )), fmt="%.2f")
    print("Done.")

def parse_args():
    from config import config
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-b", "--bathymetry", default=config["base_bathymetry"], type=Path)
    parser.add_argument("-l", "--lake_alt", default=config["lake_alt"], type=float)
    parser.add_argument("-s", "--flood_seed", default=config["flood_seed"], type=float, nargs=2)
    parser.add_argument("-d", "--dilation_radius", default=config["flood_dilation"], type=int)
    parser.add_argument("-p", "--plot", action="store_true")
    return parser.parse_args()

def main():
    args = parse_args()
    args.flood_seed = float(args.flood_seed[0]), float(args.flood_seed[1])

    x, y, Z = topo_utils.read_raster(Path(args.bathymetry).expanduser())
    x, y, Z, qinit_bbox = make_qinit(x, y, Z, args.lake_alt, args.flood_seed, args.dilation_radius)
    np.savetxt("qinit.bbox", qinit_bbox)
    save_qinit(x, y, Z)

    if args.plot:
        plt.imshow(Z, extent=(x.min(), x.max(), y.min(), y.max()))
        plt.show()


if __name__ == "__main__":
    main()
