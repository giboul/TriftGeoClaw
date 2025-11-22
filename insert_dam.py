#!/usr/bin/env python
"""
Read the `dam_shapefiles`, their altitudes and include them in `bathymetry`, 
a copy of `base_bathymetry`.
"""
from pathlib import Path
import numpy as np
from numpy.typing import ArrayLike
import topo_utils
from matplotlib.path import Path as mPath
from matplotlib import pyplot as plt


def insert_dam(dams: ArrayLike,
               dam_alts: ArrayLike,
               x: ArrayLike,
               y: ArrayLike,
               Z: ArrayLike):
    
    if not np.iterable(dam_alts):
        dam_alts = [dam_alts]

    X, Y = np.meshgrid(x, y)
    shape = Z.shape
    Z = Z.flatten()
    X = X.flatten()
    Y = Y.flatten()
    XY = np.column_stack((X, Y))

    bboxes = []

    for dam, dam_alt in zip(dams, dam_alts):
        in_dam = mPath(dam).contains_points(XY)
        Z[in_dam] = np.maximum(Z[in_dam], dam_alt)
        bboxes.append([
            X[in_dam].min(), X[in_dam].max(),
            Y[in_dam].min(), Y[in_dam].max()
        ])

    return Z.reshape(shape), bboxes


def parse_args():
    from config import config
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--dams", "-d", type=Path, default=config["dams"])
    parser.add_argument("--dam_alts", "-a", type=float, default=config["dam_alts"], nargs="?")
    parser.add_argument("--base_bathymetry", "-B", type=Path, default=config["base_bathymetry"])
    parser.add_argument("--bathymetry", "-b", type=Path, default=config["bathymetry"])
    parser.add_argument("--no_plot", "-s", action="store_true")
    parser.add_argument("--bbox", "-x", nargs=4, type=float, default=config["lower"]+config["upper"])
    parser.add_argument("--cell_size", "-c", type=float, default=config["cell_size"])
    return parser.parse_args()


def main():

    args = parse_args()

    dams = topo_utils.read_poly(Path(args.dams).expanduser())
    x, y, Z = topo_utils.read_raster(args.base_bathymetry, args.bbox, args.cell_size)

    Z, bboxes = insert_dam(dams, args.dam_alts, x, y, Z)

    bbox = x.min(), y.min(), x.max(), y.max()
    topo_utils.write_raster(Path(args.bathymetry).expanduser(), Z, bbox)
    for i, bbox in enumerate(bboxes):
        np.savetxt(f"dam{i}.bbox", bbox)

    if args.no_plot is False:
        plt.imshow(Z, extent=(x.min(), x.max(), y.min(), y.max()))
        # plt.fill(dams[0][:, 0], dams[0][:, 1])
        plt.show()


if __name__ == "__main__":
    main()
