#!/usr/bin/env python
from pathlib import Path
from argparse import ArgumentParser
from shutil import rmtree
import params
import numpy as np
from matplotlib import pyplot as plt
from tifffile import TiffFile
from scipy.interpolate import RegularGridInterpolator


def write_topo():
    
    # Temporary directory
    ifile = Path("..") / "swissALTI3D_merged.tif"
    tempdir = Path("_temp")
    tempdir.mkdir(exist_ok=True)

    print(f"\tINFO: Opening {ifile}... ")
    # With just tifffile
    with TiffFile(ifile) as tif:
        Z = tif.asarray()
        nx = tif.pages[0].tags["ImageWidth"].value
        ny = tif.pages[0].tags["ImageLength"].value
        x_res, y_res, z_res = tif.pages[0].tags["ModelPixelScaleTag"].value
        xmin = tif.pages[0].tags["ModelTiepointTag"].value[3]
        ymax = tif.pages[0].tags["ModelTiepointTag"].value[4]
    x = xmin + (np.arange(nx) + 0.5)*x_res
    y = ymax - (np.arange(ny) + 0.5)*y_res
    # Add some room to avoid interpolation error
    xmin, xmax, ymin, ymax = expand_bounds(**params.bounds)
    print(f"Expanding bounds to {xmin, ymin, xmax, ymax = }")
    xmask = (xmin <= x) & (x <= xmax)
    ymask = (ymin <= y) & (y <= ymax)
    x = x[xmask]
    y = y[ymask]
    Z = Z[ymask, :][:, xmask]
    plt.imshow(Z, extent=(x.min(), x.max(), y.min(), y.max()))
    X, Y = np.meshgrid(x, y)
    print(f"\tINFO: Downsampling and cropping {ifile} to {tempdir / 'bathy.tiff'}")
    interp = RegularGridInterpolator((x, y), Z.T)
    print(GridInterpolator(x, y, Z, x, y))
    x = np.interp(np.arange(xmin, xmax, step=params.resolution), x, x)
    y = np.interp(np.arange(ymin, ymax, step=params.resolution), y[::-1], y[::-1])[::-1]
    X, Y = np.meshgrid(x, y)
    Z = interp((X, Y))
 
    print(f"\tINFO: Adding dam... ")
    insert_dam(X, Y, Z)

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
    yt = yt[::-1]
    y = y[::-1]
    x = np.maximum(xt[0], np.minimum(xt[-1], x))
    y = np.maximum(yt[0], np.minimum(yt[-1], y))
    ix = np.searchsorted(xt, x)
    iy = np.searchsorted(yt, y)
    print(f"{ix = }")
    print(f"{iy = }")
    Zt = np.hstack((Zt[:, :1], Zt, Zt[:, -1:]))
    Zt = np.vstack((Zt[:1, :], Zt, Zt[-1:, :]))
    print(f"{Zt.shape = }")
    xslope = (Zt[:, :-2] - Zt[:, 2:])[:, ix]
    yslope = (Zt[:-2, :] - Zt[2:, :])[iy, :]
    dx = x - xt[ix-1]
    dy = y - yt[iy-1]
    print(f"{Zt.shape = }")
    np.take_along_axis(Zt, iy, ix)
    Z = Zt[1:-1, 1:-1][iy, :][:, ix] + xslope*dx + yslope*dy
    return np.take_along_axis(Z)


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

