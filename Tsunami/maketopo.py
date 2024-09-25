#!/usr/bin/env python
from pathlib import Path
from shutil import rmtree
from osgeo.gdal import UseExceptions, Open, Translate, Warp
import numpy as np
from matplotlib import pyplot as plt
from argparse import ArgumentParser
import params
from skimage.morphology import flood
# from clawpack.geoclaw.marching_front import select_by_flooding


def main():
    UseExceptions()
    
    # Temporary directory
    ifile = Path("..") / "swissALTI3D_merged.tif"
    tempdir = Path("_temp")
    tempdir.mkdir(exist_ok=True)

    if True:
        ulx, uly, lrx, lry = params.xmin, params.ymax, params.xmax, params.ymin
    else:
        ulx = params.xmin - 0.5*(params.xmax - params.xmin)
        uly = params.ymax + 0.5*(params.ymax - params.ymin)
        lrx = params.xmax + 0.5*(params.xmax - params.xmin)
        lry = params.ymin - 0.5*(params.ymax - params.ymin)

    print(f"\tINFO: Opening {ifile}... ")
    data = Open(str(ifile))

    print(f"\tINFO: Downsampling {ifile} to {tempdir / 'b1.tiff'}")
    data = Warp(str(tempdir / 'b1.tif'), data, xRes=params.xRes, yRes=params.yRes)
    
    print(f"\tINFO: Cropping {tempdir / 'b1.tif'} to {tempdir / 'b2.tif'}... ")
    data = Translate(str(tempdir / "b2.tif"), data, projWin=(ulx, uly, lrx, lry))

    print(f"\tINFO: Converting {tempdir / 'b2.tif'} to bathymetry.xyz... ")
    Translate(str(tempdir / "bathymetry.xyz"), data, format="xyz")
    
    print(f"\tINFO: Drawing dam... ")
    x, y, z = np.loadtxt(tempdir / "bathymetry.xyz").T
    dam_y1 = params.dam_upstream(x, y)
    dam_y2 = params.dam_downstream(x, y) 
    z[(dam_y1 <= y) & (y <= dam_y2) & (z < params.dam_z)] = params.dam_z
    print(f"\tINFO: Saving bathy_with_dam.xyz... ")
    np.savetxt("bathy_with_dam.xyz", np.vstack((x, y, z)).T)
    print(f"\t\tFile size is {Path('bathy_with_dam.xyz').stat().st_size:.2g} bytes")

    print("\tINFO: Writing qinit.xyz... ")
    ny = np.unique(y).size
    nx = y.size // ny
    z_lake = z.reshape(ny, nx).copy()
    seed_x = 2.6703e6
    seed_y = 1.1715e6
    seed = ((x-seed_x)**2 + (y-seed_y)**2).argmin()
    seed = seed//nx, seed%nx
    print(f"\t\tFlooding around {seed_x} {seed_y}...")
    flooded = flood(z_lake < params.lake_level, seed)
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
    main()

