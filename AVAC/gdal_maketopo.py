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
    np.savetxt(f"bathymetry.xyz", np.vstack((x, y, z)).T)

    # Adding dam
    # dam = ((x >= 2.6700500) & (x <= 2.6705) &
    #        (y >= 1.1717714) & (y <= 1.1719301))
    # z[dam] += 115  # m (hauteur du barrage)

    # print(f"\tINFO: Converting {tempdir / 'b2.tif'} to bathymetry.asc")
    # Translate(f"bathymetry.asc", data, format="AAIGrid")  # Not in the right format

    rmtree(tempdir)


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
