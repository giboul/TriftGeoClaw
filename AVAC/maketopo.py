#!/usr/bin/env python
"""
Script for generating the topography file and the extents.

Currently
---------
    `gdaltranslate` crops, coarsens and the converts the file to `bathymetry.xyz`
    `plot` helps to analyse both created files
    `make_qinit` as well as `qinit` are in progress
"""
from pathlib import Path
from shutil import rmtree
from osgeo.gdal import Open, Translate, Warp

import numpy as np
from clawpack.geoclaw.topotools import Topography

from matplotlib import pyplot as plt
from AddSetrun import (
    xmin, xmax,
    ymin, ymax,
    sea_level
)


def gdaltranslate() -> None:

    xRes = yRes = 10
    # Temporary directory
    ifile = Path("..")/"swissALTI3D_merged.tif"
    tempdir = Path("_temp")
    tempdir.mkdir(exist_ok=True)

    ulx, uly, lrx, lry = xmin, ymax, xmax, ymin

    data = Open(str(ifile))
    print(f"\tINFO: Cropping {ifile} to {tempdir / 'b1.tif'}")
    data = Translate(str(tempdir / "b1.tif"), data, projWin=(ulx, uly, lrx, lry))
    print(f"\tINFO: Downsampling {tempdir / 'b1.tif'} to {xRes=} {yRes=}")
    data = Warp(str(tempdir / 'b2.tif'), data, xRes=xRes, yRes=yRes)

    print(f"\tINFO: Converting {tempdir / 'b2.tif'} to bathymetry.xyz")
    Translate(f"bathymetry.xyz", data, format="xyz")

    print(f"\tINFO: Rescaling bathymetry.xyz to computable values")
    x, y, z = np.loadtxt("bathymetry.xyz").T
    print(f"\tINFO: Stripping lake level to altitude 0 (bathymetry.xyz)")
    z = z - sea_level
    np.savetxt(f"bathymetry.xyz", np.vstack((x, y, z)).T)

    # Adding dam
    # dam = ((x >= 2.6700500) & (x <= 2.6705) &
    #        (y >= 1.1717714) & (y <= 1.1719301))
    # z[dam] += 115  # m (hauteur du barrage)

    # print(f"\tINFO: Converting {tempdir / 'b2.tif'} to bathymetry.asc")
    # Translate(f"bathymetry.asc", data, format="AAIGrid")  # Not in the right format

    rmtree(tempdir)


def plot(asc=False, h0=False) -> None:

    print(f"\tINFO: plotting bathymetry and approximate computational domain")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    if asc is True:
        nx, ny, x0, y0, cs = np.loadtxt("bathymetry.asc", max_rows=5, usecols=1)
        nx, ny = int(nx), int(ny)
        x = np.linspace(x0, x0 + nx*cs, num=nx)
        y = np.linspace(y0, y0 + ny*cs, num=ny)
        z = np.loadtxt("bathymetry.asc", skiprows=5)
        x, y = np.meshgrid(x, y)
    else:
        x, y, z = np.loadtxt("bathymetry.xyz").T

    ax.plot_trisurf(x[::20], y[::20], z[::20], color='gray', alpha=0.7, label="Bathymetry")

    if h0 is True:
        xw, yw, zw = np.loadtxt("qinit.xyz").T
        above = zw >= z
        # above = 0 >= z
        if above.sum() >= 3:
            ax.plot_trisurf(x[above], y[above], np.full(above.sum(), 0), label="q0")

    # Plot approximate computational domain limits
    walls_kwargs = dict(facecolor='red', alpha=.2)
    Xc, Zc = np.meshgrid((xmin, xmax), (z.min(), z.max()))
    Yc = np.full_like(Xc, ymin)
    ax.plot_surface(Xc, Yc, Zc, **walls_kwargs, label="Computational bounds")

    Xc, Zc = np.meshgrid((xmin, xmax), (z.min(), z.max()))
    Yc = np.full_like(Xc, ymax)
    ax.plot_surface(Xc, Yc, Zc, **walls_kwargs)

    Yc, Zc = np.meshgrid((ymin, ymax), (z.min(), z.max()))
    Xc = np.full_like(Xc, xmin)
    ax.plot_surface(Xc, Yc, Zc, **walls_kwargs)

    Yc, Zc = np.meshgrid((ymin, ymax), (z.min(), z.max()))
    Xc = np.full_like(Xc, xmax)
    ax.plot_surface(Xc, Yc, Zc, **walls_kwargs)

    ax.set_aspect("equal")
    ax.legend(loc="upper left")
    if True:
        file = "Topography.png"
        print(f"\tINFO: Saving {file}...", end=" ")
        fig.savefig(file, bbox_inches='tight')
        print("Done.")
    else:
        plt.show()


def make_qinit():
    """
    Create qinit data file
    """
    x, y, z = np.loadtxt("bathymetry.xyz").T
    nxpoints = np.unique(x).size
    nypoints = np.unique(y).size
    xlower = x.min()
    xupper = x.max()
    ylower = y.min()
    yupper = y.max()
    outfile = "qinit.xyz"

    topography = Topography(topo_func=qinit)
    topography.x = np.linspace(xlower,xupper,nxpoints)
    topography.y = np.linspace(ylower,yupper,nypoints)
    topography.write(outfile)


def qinit(x, y):
    """"""
    q = np.zeros_like(x+y, dtype=np.float16)
    # q[(np.abs(x-cxstart) <= 1e-3) & (np.abs(y-1.1714) <= 1e-3)] = 100
    xm = 2.6690e6  # (xstart+xstop)/2
    ym = 1.1704e6  # 1171400.-ystart  # (ystart+ystop)/2
    # q += sea_level
    domain = (x-xm)**2+(y-ym)**2 <= (2e2)**2
    q[domain] = 10
    # print(q.min(), q.max())
    return q


if __name__ == "__main__":
    print(f"INFO: running {__file__}")
    gdaltranslate()
    make_qinit()
    plot(h0=True)
