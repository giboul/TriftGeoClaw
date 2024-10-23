#!/usr/bin/env python
from pathlib import Path
from argparse import ArgumentParser
from shutil import rmtree
import AddSetrun
import numpy as np
from matplotlib import pyplot as plt
from skimage.morphology import flood, isotropic_dilation
from skimage.measure import find_contours
from tifffile import TiffFile
from scipy.interpolate import RegularGridInterpolator as RGI
# from osgeo.gdal import UseExceptions, Open, Translate, Warp
# from clawpack.geoclaw.marching_front import select_by_flooding
params = AddSetrun

def write_topo(plot=False):

    # Temporary directory
    ifile = Path("..") / "swissALTI3D_merged.tif"
    tempdir = Path("_temp")
    tempdir.mkdir(exist_ok=True)

    print(f"\tINFO: Opening {ifile}... ")
    # With just tifffile
    with TiffFile(ifile) as tif:
        Ztif = tif.asarray()
        nx = tif.pages[0].tags["ImageWidth"].value
        ny = tif.pages[0].tags["ImageLength"].value
        x_res, y_res, z_res = tif.pages[0].tags["ModelPixelScaleTag"].value
        xmin = tif.pages[0].tags["ModelTiepointTag"].value[3]
        ymax = tif.pages[0].tags["ModelTiepointTag"].value[4]
    if plot:
        plt.title("Original topography")
        plt.imshow(Ztif)
        plt.show()
    xtif = xmin + (np.arange(nx) + 0.5)*x_res
    ytif = ymax - (np.arange(ny) + 0.5)*y_res
    # Add some room to avoid interpolation error
    xmin, xmax, ymin, ymax = expand_bounds(**params.bounds)

    print(f"\tCropping to {xmin, ymin, xmax, ymax = }")
    xmask = (xmin <= xtif) & (xtif <= xmax)
    ymask = (ymin <= ytif) & (ytif <= ymax)
    xtif = xtif[xmask]
    ytif = ytif[ymask][::-1]
    Ztif = Ztif[ymask, :][:, xmask]

    print(f"\tDownscaling to resolution = {params.resolution}")
    x = np.arange(xmin, xmax+params.resolution, step=params.resolution)
    y = np.arange(ymin, ymax+params.resolution, step=params.resolution)
    # X, Y = np.meshgrid(x, y)
    # Z = RGI((ytif, xtif), Ztif, fill_value=None, bounds_error=False)
    # Z = Z(np.vstack((Y.flatten(), X.flatten())).T).reshape(y.size, x.size)
    Z = grid_interp(xtif, ytif, Ztif, x, y)

    print(f"\tINFO: Saving bathy_with_dam.asc... ")
    asc_header = "\n".join((
        f"{x.size} ncols",
        f"{y.size} nrows",
        f"{x.min()} xllcenter",
        f"{y.min()} yllcenter",
        f"{params.resolution} cellsize",
        f"{999999} nodata_value"
    ))
    np.savetxt("bathy_with_dam.asc", Z.flatten(), header=asc_header, comments="")
    print(f"\tINFO: File size is {Path('bathy_with_dam.asc').stat().st_size:.2g} bytes.")

    rmtree(tempdir)

    if plot:
        extent = (xmin, xmax, ymin, ymax) 
        # h = params.lake_alt - Z
        # h[h <= 0] = float("nan")
        plt.title("Processed topography")
        plt.imshow(Z, extent=extent)
        # plt.imshow(h, cmap="Blues", extent=extent)
        plt.show()


def expand_bounds(xmin, xmax, ymin, ymax, margin=0.25):
    _xmin = xmin - margin*(xmax - xmin)
    _ymin = ymin - margin*(ymax - ymin)
    _xmax = xmax + margin*(xmax - xmin)
    _ymax = ymax + margin*(ymax - ymin)
    return _xmin, _xmax, _ymin, _ymax


def grid_interp(xt, yt, Zt, x, y):
    x = np.clip(x, xt.min(), xt.max())
    y = np.clip(y, yt.min(), yt.max())
    ix = np.clip(np.searchsorted(xt, x)-1, 0, xt.size-1)
    iy = np.clip(np.searchsorted(yt, y)-1, 0, yt.size-1)
    Z11 = Zt[iy, :][:, ix]
    Z21 = Zt[iy, :][:, ix+1]
    Z12 = Zt[iy+1, :][:, ix]
    Z22 = Zt[iy+1, :][:, ix+1]
    x1 = xt[ix]
    x2 = xt[ix+1]
    y1 = yt[iy]
    y2 = yt[iy+1]
    Z = ((
        + (x2-x)*((y2-y)*Z11.T).T
        + (x2-x)*((y-y1)*Z12.T).T
        + (x-x1)*((y2-y)*Z21.T).T
        + (x-x1)*((y-y1)*Z22.T).T
    ).T/(y2-y1)).T/(x2-x1)
    return Z


def write_topo_old():
    UseExceptions()
    
    # Temporary directory
    ifile = Path("..") / "swissALTI3D_merged.tif"
    tempdir = Path("_temp")
    tempdir.mkdir(exist_ok=True)

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
    z[(dam_y1 <= y) & (y <= dam_y2) & (z < params.dam_level)] = params.dam_level
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

    print(f"\tINFO: Writing lake countour...")
    seed_x, seed_y = params.flood_seed
    seed = ((x-seed_x)**2 + (y-seed_y)**2).argmin()
    seed = seed//nx, seed%nx
    print(f"\t\tFlooding around {seed_x} {seed_y}...")
    flooded = fill_lake(z.reshape(ny, nx).copy(), seed, params.lake_level, 10/params.resolution)
    print(f"\t\tFinding flood bounding box...")
    yc, xc = find_contours(flooded, 0.5)[0].T
    xc = xmin + xc/nx * (xmax - xmin)
    yc = ymin + yc/ny * (ymax - ymin)
    yc = ymin + (ymax-yc)
    with open("lake_bounds.py", "w") as file:
        file.write("\n".join((
            f"xmin = {xc.min()}",
            f"xmax = {xc.max()}",
            f"ymin = {yc.min()}",
            f"ymax = {yc.max()}"
        )))

    rmtree(tempdir)

    parser = ArgumentParser()
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    if args.plot or True:
        extent =(xmin, xmax, ymin, ymax) 
        plt.imshow(z.reshape(ny, nx), extent=extent)
        xcmin, xcmax = xc.min(), xc.max()
        ycmin, ycmax = yc.min(), yc.max()
        plt.plot((xcmin, xcmax, xcmax, xcmin, xcmin),
                 (ycmin, ycmin, ycmax, ycmax, ycmin))
        plt.plot(xc, yc, label="lake contour", c="g")
        plt.scatter(seed_x, seed_y, label="flood seed", c='k')
        plt.legend()
        plt.show()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    write_topo(plot=args.plot)

