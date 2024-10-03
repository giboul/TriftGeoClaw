#!/usr/bin/env python
from pathlib import Path
from argparse import ArgumentParser
from shutil import rmtree
import params
import numpy as np
from matplotlib import pyplot as plt
from skimage.morphology import flood, isotropic_dilation
from skimage.measure import find_contours
from osgeo.gdal import UseExceptions, Open, Translate, Warp
# from clawpack.geoclaw.marching_front import select_by_flooding


def write_topo():
    UseExceptions()
    
    # Temporary directory
    ifile = Path("..") / "swissALTI3D_merged.tif"
    tempdir = Path("_temp")
    tempdir.mkdir(exist_ok=True)

    bds = params.bounds  # Add some room for error
    xmin = bds["xmin"] - (bds["xmax"] - bds["xmin"])/2
    ymin = bds["ymin"] - (bds["ymax"] - bds["ymin"])/2
    xmax = bds["xmax"] + (bds["xmax"] - bds["xmin"])/2
    ymax = bds["ymax"] + (bds["ymax"] - bds["ymin"])/2

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
    dam_y1 = dam_upstream(x, y)
    dam_y2 = dam_downstream(x, y) 
    z[(dam_y1 <= y) & (y <= dam_y2) & (z < params.dam_alt)] = params.dam_alt
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
    np.savetxt("bathy_with_dam.asc", z, header=asc_header)
    print(f"\t\tFile size is {Path('bathy_with_dam.asc').stat().st_size:.2g} bytes.")

    print("\tINFO: Writing qinit.xyz... ")
    z_lake = z.reshape(ny, nx).copy()
    seed_x, seed_y = params.flood_seed
    seed = ((x-seed_x)**2 + (y-seed_y)**2).argmin()
    seed = seed//nx, seed%nx
    print(f"\t\tFlooding around {seed_x} {seed_y}...")
    flooded = fill_lake(z_lake, seed, params.lake_alt, 10/params.resolution)
    z_lake = z_lake.flatten()
    np.savetxt("qinit.xyz", np.vstack((x, y, z_lake)).T)
    print(f"\t\tFile size is {Path('qinit.xyz').stat().st_size:.2g} bytes.")
    print(f"\tINFO: Writing lake countour...")
    yc, xc = find_contours(flooded, 0.5)[0].T
    xc = xmin + xc/nx * (xmax - xmin)
    yc = ymin + yc/ny * (ymax - ymin)
    yc = ymin + (ymax-yc)
    np.savetxt("lake_contour.xy", np.vstack((xc, yc)).T)
    print(f"\t\tFile size is {Path('lake_contour.xy').stat().st_size:.2g} bytes.")

    rmtree(tempdir)

    parser = ArgumentParser()
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    if args.plot or True:
        extent =(xmin, xmax, ymin, ymax) 
        h = z_lake - z
        h[h <= 0] = float("nan")
        plt.imshow(z.reshape(ny, nx), extent=extent, cmap="inferno")
        plt.imshow(h.reshape(ny, nx), cmap="Blues", extent=extent)
        plt.plot(xc, yc, label="lake contour", c="g")
        plt.scatter(seed_x, seed_y, label="flood seed", c='k')
        plt.legend()
        plt.show()


def fill_lake(topo, seed, max_level=0, dilation_radius=0):
    flooded = flood(topo < max_level, seed)
    isotropic_dilation(flooded, dilation_radius, flooded)
    topo[flooded] = max_level
    return flooded


def dam_upstream(x, y, l=2670561, offset=0):
    bds = params.bounds
    yd = bds["ymax"] - 0.3*(x-bds["xmin"]) - 50000/(x-l+offset) - 300 + offset
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

