#!/usr/bin/env python
from pathlib import Path
from shutil import rmtree
from osgeo.gdal import UseExceptions, Open, Translate, Warp
import numpy as np
from matplotlib import pyplot as plt
from argparse import ArgumentParser
import params
# from clawpack.geoclaw.marching_front import select_by_flooding


def main():
    UseExceptions()
    
    # Temporary directory
    ifile = Path("..") / "swissALTI3D_merged.tif"
    tempdir = Path("_temp")
    tempdir.mkdir(exist_ok=True)
    
    ulx, uly, lrx, lry = params.xmin, params.ymax, params.xmax, params.ymin
    
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

    mask = (dam_y1 <= y) & (y <= dam_y2) & (z < params.dam_z)
    z[mask] = params.dam_z

    print(f"\tINFO: Saving bathy_with_dam.xyz... ")
    np.savetxt("bathy_with_dam.xyz", np.vstack((x, y, z)).T)
    print(f"\t\tFile size is {Path('bathy_with_dam.xyz').stat().st_size:.2g} bytes")

    print("\tINFO: Writing qinit.xyz... ")
    z_lake = np.full_like(z, params.lake_level)
    no_alt = 0
    z_lake[y > dam_y2] = no_alt
    z_lake[
      (x < 2669800) | (2670900 < x) | (y < 1170300) | (1172000 < y)
    ] = no_alt
    np.savetxt("qinit.xyz", np.vstack((x, y, z_lake)).T)
    print(f"\t\tFile size is {Path('qinit.xyz').stat().st_size:.2g} bytes")

    rmtree(tempdir)
    
    parser = ArgumentParser()
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    if args.plot:
        ny = np.unique(y).size
        nx = y.size // ny
        extent = (ulx, lrx, lry, uly)
        print(f"{extent = }")
        plt.imshow(z.reshape(ny, nx), extent=extent, cmap="inferno")
        h = z_lake - z
        h[h <= 0] = float("nan")
        plt.imshow(h.reshape(ny, nx), cmap="Blues")
        plt.gca().set_aspect("equal")
        plt.show()


if __name__ == "__main__":
    main()

