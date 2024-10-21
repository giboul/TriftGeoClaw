#!/usr/bin/env python
from pathlib import Path
from argparse import ArgumentParser
from shutil import rmtree
import params
import numpy as np
from matplotlib import pyplot as plt
from osgeo.gdal import UseExceptions, Open, Translate, Warp
# from clawpack.geoclaw.marching_front import select_by_flooding


def write_topo():
    UseExceptions()
    
    # Temporary directory
    ifile = Path("..") / "swissALTI3D_merged.tif"
    tempdir = Path("_temp")
    tempdir.mkdir(exist_ok=True)

    # Add some room to avoid interpolation error
    xmin, xmax, ymin, ymax = expand_bounds(**params.bounds)
    print(f"Expanding bounds to {xmin, ymin, xmax, ymax = }")

    print(f"\tINFO: Opening {ifile}... ")
    data = Open(str(ifile))

    print(f"\tINFO: Downsampling and cropping {ifile} to {tempdir / 'bathy.tiff'}")
    data = Warp(str(tempdir / 'bathy.tif'), data,
                xRes=params.resolution,
                yRes=params.resolution,
                outputBounds=(xmin, ymin, xmax, ymax))
 
    print(f"\tINFO: Converting {tempdir / 'bathy.tif'} to bathy.xyz... ")
    Translate(str(tempdir / "bathy.xyz"), data, format="xyz")

    print(f"\tINFO: Drawing dam... ")
    x, y, z = np.loadtxt(tempdir / "bathy.xyz").T
    insert_dam(x, y, z)

    print(f"\tINFO: Saving bathy_with_dam.asc... ")
    ny = np.unique(y).size
    nx = y.size // ny
    asc_header = "\n".join((
        f"{nx} ncols",
        f"{ny} nrows",
        f"{x.min()} xllcenter",
        f"{y.min()} yllcenter",
        f"{params.resolution} cellsize",
        f"{999999} nodata_value"
    ))
    np.savetxt("bathy_with_dam.asc", z, header=asc_header, comments="")
    print(f"\tINFO: File size is {Path("bathy_with_dam.asc").stat().st_size:.2g} bytes.")

    rmtree(tempdir)

    parser = ArgumentParser()
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    if args.plot:
        extent = (xmin, xmax, ymin, ymax) 
        h = params.lake_alt - z
        h[h <= 0] = float("nan")
        plt.imshow(z.reshape(ny, nx), extent=extent)
        plt.imshow(h.reshape(ny, nx), cmap="Blues", extent=extent)
        plt.show()


def expand_bounds(xmin, xmax, ymin, ymax, margin=0.5):
    _xmin = xmin - margin*(xmax - xmin)
    _ymin = ymin - margin*(ymax - ymin)
    _xmax = xmax + margin*(xmax - xmin)
    _ymax = ymax + margin*(ymax - ymin)
    return _xmin, _xmax, _ymin, _ymax


def insert_dam(x, y, z):
    dam_y1 = dam_upstream(x, y)
    dam_y2 = dam_downstream(x, y) 
    z[(dam_y1 <= y) & (y <= dam_y2) & (z < params.dam_alt)] = params.dam_alt


def dam_upstream(x, y, offset=0, y0=1171960, x0=2669850, x1=2670561):
    yd = y0 - 0.3*(x+offset*0-x0) - 50000/(x+offset*0-x1) - 300 + offset
    yd[(x < x0) | (x > x1)] = float("inf")
    return yd


def dam_downstream(x, y, thk=30):
    d = dam_upstream(x, y)
    u = dam_upstream(x+thk, y, offset=30)
    return np.maximum(u, d)


if __name__ == "__main__":
    write_topo()

