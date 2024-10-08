#!/usr/bin/env python
from pathlib import Path
from argparse import ArgumentParser
from shutil import rmtree
import params
import numpy as np
from matplotlib import pyplot as plt
from tifffile import TiffFile


def write_topo():
    
    # Temporary directory
    ifile = Path("..") / "swissALTI3D_merged.tif"
    ifile = "mock_topo.tif"
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
    xtif = xmin + (np.arange(nx) + 0.5)*x_res
    ytif = ymax - (np.arange(ny) + 0.5)*y_res
    # Add some room to avoid interpolation error
    xmax, ymin = xtif.max(), ytif.min()  # xmin, xmax, ymin, ymax = expand_bounds(**params.bounds)

    print(f"\tCropping to {xmin, ymin, xmax, ymax = }")
    xmask = (xmin <= xtif) & (xtif <= xmax)
    ymask = (ymin <= ytif) & (ytif <= ymax)
    xtif = xtif[xmask]
    ytif = ytif[ymask][::-1]
    Ztif = Ztif[ymask, :][:, xmask]
    
    print(f"\tDownscaling to resolution = {params.resolution}")
    x = np.arange(xmin, xmax+params.resolution, step=params.resolution)
    y = np.arange(ymin, ymax+params.resolution, step=params.resolution)
    X, Y = np.meshgrid(x, y)
    Z = GridInterpolator(xtif, ytif, Ztif, x, y)
 
    print(f"\tINFO: Adding dam... ")
    Z = insert_dam(X, Y, Z)

    print(f"\tINFO: Saving bathy_with_dam.asc... ")
    asc_header = "\n".join((
        f"{nx} ncols",
        f"{ny} nrows",
        f"{x.min()} xllcenter",
        f"{y.min()} yllcenter",
        f"{params.resolution} cellsize",
        f"{999999} nodata_value"
    ))
    np.savetxt("bathy_with_dam.asc", Z.flatten(), header=asc_header, comments="")
    print(f"\tINFO: File size is {Path("bathy_with_dam.asc").stat().st_size:.2g} bytes.")

    rmtree(tempdir)

    parser = ArgumentParser()
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    if args.plot or True:
        extent = (xmin, xmax, ymin, ymax) 
        h = params.lake_alt - Z
        h[h <= 0] = float("nan")
        plt.imshow(Z, extent=extent, cmap="inferno", vmin=params.lake_alt)
        plt.imshow(h, cmap="Blues", extent=extent)
        plt.show()


def GridInterpolator(xt, yt, Zt, x, y):
    # Force input array inside bounds
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharex=True, sharey=True)
    ax1.imshow(Zt, extent=(xt.min(), xt.max(), yt.min(), yt.max()), cmap="inferno")
    x = np.maximum(xt[0], np.minimum(xt[-1], x))
    y = np.maximum(yt[0], np.minimum(yt[-1], y))
    # Find nearest indices
    ix = np.searchsorted(xt, x)
    iy = np.searchsorted(yt, y)
    # Compute slopes
    Zt = np.hstack((Zt, Zt[:, -1:]))
    Zt = np.vstack((Zt, Zt[-1:, :]))
    dx = xt[1:] - xt[:-1]
    dy = yt[1:] - yt[:-1]
    dx = np.hstack((dx, dx[-1]))
    dy = np.hstack((dy, dy[-1]))
    xslope = ((Zt[:, 1:] - Zt[:, :-1])/dx)[:-1, :]
    yslope = ((Zt[1:, :] - Zt[:-1, :]).T/dy).T[:, :-1]
    xslope = xslope[iy, :][:, ix]
    yslope = yslope[iy, :][:, ix]
    # Remove the extension
    Z = Zt[:-1, :-1][iy, :][:, ix]
    ax2.imshow(xslope, extent=(x.min(), x.max(), y.min(), y.max()), cmap="inferno")
    # Add the increments
    Z = Z + (x - xt[ix])*xslope + ((y - yt[iy])*yslope.T).T
    ax3.imshow((xt[ix]-x)*xslope, extent=(x.min(), x.max(), y.min(), y.max()), cmap="inferno")
    plt.show()
    exit()


def expand_bounds(xmin, xmax, ymin, ymax, margin=0.5):
    _xmin = xmin - margin*(xmax - xmin)
    _ymin = ymin - margin*(ymax - ymin)
    _xmax = xmax + margin*(xmax - xmin)
    _ymax = ymax + margin*(ymax - ymin)
    return _xmin, _xmax, _ymin, _ymax


def insert_dam(x, y, z):
    y = y[::-1]
    dam_y1 = dam_upstream(x, y)
    dam_y2 = dam_downstream(x, y) 
    z[(dam_y1 <= y) & (y <= dam_y2) & (z < params.dam_alt)] = params.dam_alt
    return z


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

