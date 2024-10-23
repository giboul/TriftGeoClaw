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
    xtif = xmin + (np.arange(nx) + 0.5)*x_res
    ytif = ymax - (np.arange(ny) + 0.5)*y_res
    if plot:
        extent = (xmin, xtif.max(), ytif.min(), ymax) 
        fig, (ax1, ax2) = plt.subplots(ncols=2)
        ax1.set_title("Original topography")
        ax1.imshow(Ztif, extent=extent)
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
    Z = grid_interp(xtif, ytif, Ztif, x, y)
    Z = insert_dam(*np.meshgrid(x, y), Z)
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

    y = y[::-1]
    seed = (np.abs(x-params.flood_seed[0]).argmin(),
            np.abs(y-params.flood_seed[1]).argmin())
    flooded = flood(Z <= params.lake_alt, seed[::-1])
    xc, yc = find_contours(flooded.T, 0.5)[0].T
    xc, yc = find_contours(flooded.T, 0.5)[0].T
    xc = xmin + xc/x.size * (xmax - xmin)
    yc = ymin + yc/y.size * (ymax - ymin)
    yc = ymin + (ymax-yc)  # Reverse y
    contour_coords = np.vstack((xc, yc)).T
    np.savetxt("contour.xy", contour_coords)
    dilated = isotropic_dilation(flooded, 20/params.resolution)
    xc, yc = find_contours(dilated.T, 0.5)[0].T
    xc, yc = find_contours(dilated.T, 0.5)[0].T
    xc = xmin + xc/x.size * (xmax - xmin)
    yc = ymin + yc/y.size * (ymax - ymin)
    yc = ymin + (ymax-yc)  # Reverse y
    dilated_contour_coords = np.vstack((xc, yc)).T
    np.savetxt("contour_dilated.xy", dilated_contour_coords)

    if plot:
        extent = (xmin, xmax, ymin, ymax) 
        ax2.set_title("Processed topography")
        ax2.imshow(Z, extent=extent)
        ax2.plot(*contour_coords.T, '-')
        ax2.plot(*dilated_contour_coords.T, '-')
        plt.scatter(x[seed[0]], y[seed[1]], c="g")
        plt.show()


def expand_bounds(xmin, xmax, ymin, ymax, margin=5*params.resolution):
    _xmin = xmin - margin
    _ymin = ymin - margin
    _xmax = xmax + margin
    _ymax = ymax + margin
    return _xmin, _xmax, _ymin, _ymax


def grid_interp(xt, yt, Zt, x, y):
    x = np.clip(x, xt.min(), xt.max())
    y = np.clip(y, yt.min(), yt.max())
    ix = np.clip(np.searchsorted(xt, x)-1, 0, xt.size-2)
    iy = np.clip(np.searchsorted(yt, y)-1, 0, yt.size-2)
    iyz = yt.size-2-iy
    Z11 = Zt[iyz, :][:, ix]
    Z21 = Zt[iyz, :][:, ix+1]
    Z12 = Zt[iyz+1, :][:, ix]
    Z22 = Zt[iyz+1, :][:, ix+1]
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
    return Z[::-1, :]


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


def insert_dam(x, y, z):
    dam_y1 = dam_upstream(x)
    dam_y2 = dam_downstream(x)
    y = y[::-1]
    z[(dam_y1 <= y) & (y <= dam_y2) & (z < params.dam_alt)] = params.dam_alt
    return z


def dam_upstream(x, offset=0, y0=1171960, x0=2669850, x1=2670561):
    yd = y0 - 0.3*(x+offset*0-x0) - 50000/(x+offset*0-x1) - 300 + offset
    yd[(x < x0) | (x > x1)] = float("inf")
    return yd

def dam_downstream(x, thk=30):
    d = dam_upstream(x)
    u = dam_upstream(x+thk, offset=30)
    return np.maximum(u, d)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    write_topo(plot=args.plot)

