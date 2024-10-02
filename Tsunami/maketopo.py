#!/usr/bin/env python
from pathlib import Path
from argparse import ArgumentParser
from shutil import rmtree
import numpy as np
from matplotlib import pyplot as plt
from skimage.morphology import flood
from yaml import safe_load
from osgeo.gdal import UseExceptions, Open, Translate, Warp
# from clawpack.geoclaw.marching_front import select_by_flooding


with open("settings.yaml") as file:
    settings = safe_load(file)


def write_topo():
    UseExceptions()
    
    # Temporary directory
    ifile = Path("..") / "swissALTI3D_merged.tif"
    tempdir = Path("_temp")
    tempdir.mkdir(exist_ok=True)

    bds = settings["bounds"]  # Add some room for error
    xmin = bds["xmin"] - (bds["xmax"] - bds["xmin"])/2
    ymin = bds["ymax"] + (bds["ymax"] - bds["ymin"])/2
    xmax = bds["xmax"] + (bds["xmax"] - bds["xmin"])/2
    ymax = bds["ymin"] - (bds["ymax"] - bds["ymin"])/2

    print(f"\tINFO: Opening {ifile}... ")
    data = Open(str(ifile))

    print(f"\tINFO: Downsampling {ifile} to {tempdir / 'b1.tiff'}")
    data = Warp(str(tempdir / 'b1.tif'), data,
                xRes=settings["resolution"],
                yRes=settings["resolution"],
                outputBounds=(xmin, ymin, xmax, ymax))
    
    print(f"\tINFO: Cropping {tempdir / 'b1.tif'} to {tempdir / 'b2.tif'}... ")
    data = Translate(str(tempdir / "b2.tif"), data)
    print(f"\tINFO: Converting {tempdir / 'b2.tif'} to bathymetry.xyz... ")
    Translate(str(tempdir / "bathymetry.xyz"), data, format="xyz")

    print(f"\tINFO: Drawing dam... ")
    x, y, z = np.loadtxt(tempdir / "bathymetry.xyz").T
    dam_y1 = dam_upstream(x, y)
    dam_y2 = dam_downstream(x, y) 
    z[(dam_y1 <= y) & (y <= dam_y2) & (z < settigs["dam_alt"])] = settings["dam_alt"]
    print(f"\tINFO: Saving bathy_with_dam.asc... ")
    ny = np.unique(y).size
    nx = y.size // ny
    asc_header = "\n".join((
        f"{nx} ncols",
        f"{ny} nrows",
        f"{x.min()} xll",
        f"{y.min()} yll",
        f"{settings['resolution']} cellsize",
        f"{999999} nodata_value"
    ))
    np.savetxt("bathy_with_dam.asc", z, header=asc_header)
    print(f"\t\tFile size is {Path('bathy_with_dam.asc').stat().st_size:.2g} bytes")

    print("\tINFO: Writing qinit.xyz... ")
    z_lake = z.reshape(ny, nx).copy()
    seed_x = 2.67026e6
    seed_y = 1.171490e6
    seed = ((x-seed_x)**2 + (y-seed_y)**2).argmin()
    seed = seed//nx, seed%nx
    print(f"\t\tFlooding around {seed_x} {seed_y}...")
    flooded = flood(z_lake < settings["lake_level"], seed)
    z_lake[flooded] = params.lake_level
    z_lake[~flooded] = 0
    z_lake = z_lake.flatten()
    np.savetxt("qinit.xyz", np.vstack((x, y, z_lake)).T)
    print(f"\t\tFile size is {Path('qinit.xyz').stat().st_size:.2g} bytes")

    rmtree(tempdir)

    parser = ArgumentParser()
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    if args.plot or True:
        extent = (ulx, lrx, lry, uly)
        plt.imshow(z.reshape(ny, nx), extent=extent, cmap="inferno")
        h = z_lake - z
        h[h < 0] = float("nan")
        plt.imshow(h.reshape(ny, nx), cmap="Blues", extent=extent)
        plt.scatter(seed_x, seed_y, label="flood seed", c='k')
        plt.gca().set_aspect("equal")
        plt.legend()
        plt.show()


def dam_upstream(x, y, l=2670561, offset=0):
    yd = params.ymax - 0.3*(x-params.xmin) - 50000/(x-l+offset) - 300 + offset
    yd[x > l] = float("inf") 
    return yd

def dam_downstream(x, y, thk=30):
    d = dam_upstream(x, y)
    u = dam_upstream(x+thk, y, offset=30)
    return np.maximum(u, d)

def flood_mask(zimage, seed_point, max_level=0):
    plt.imshow(zimage)
    plt.gcf().show()
    below_water_level = zimage < max_level
    flooded = flood(below_water_level, seed_point)
    plt.scatter(*seed_point[::-1], c='r')
    plt.figure()
    im = plt.imshow(flooded)
    plt.colorbar(im)
    plt.show()
    return flooded


if __name__ == "__main__":
    write_topo()

