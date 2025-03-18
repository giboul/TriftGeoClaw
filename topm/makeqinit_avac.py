#!/usr/bin/env python
"""
Module for writing the qinit.xyz describing
the initial depth of the snow flow in AVAC.

Running
-------
Either call `write_qinit(avid: int=1, plot: bool=False)`
or run this script with `python makeqinit_avac.py <avalanche_id: int=1> [--plot]`
"""
from sys import getrecursionlimit, setrecursionlimit
from pathlib import Path
from yaml import safe_load
import numpy as np
from matplotlib.path import Path as mPath
from matplotlib import pyplot as plt
from skimage.morphology import isotropic_dilation
import utils


projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    TOPM = config["TOPM"]
    AVAC = config["AVAC"]


def write_qinit(avid=1, plot=False):
    """
    1. Read bathymetry.
    2. Read avalanches (geojson file) exported from from QGIS polygons shapefile.
    3. Create qinit.xyz for avalanche <avid>.
    4. If desired, plot the bathymetry and qinit.

    Note
    ----
    Only the initial depth is treated here (no initial speeds).

    Parameters
    ----------
    avid: int
        The id of the avalanche (from the geojson file)
    plot: bool = False
        Plot qinit on top of bathymetry
    """

    x, y, Z = utils.read_asc(projdir/TOPM["bathymetry"])
    geopath = projdir / TOPM["avalanches"]
    if geopath.suffix == ".geojson":
        geojson = utils.read_geojson(geopath)
    elif geopath.suffix == ".csv":
        geojson = np.loadtxt(geopath, delimiter=",").T
    else:
        geojson = np.loadtxt(geopath).T

    if avid == -1: # Launch all avalanches (for debugging purposes in TSUL)
        qinit = np.zeros(Z.shape)
        for avid in np.unique(geojson[0]):
            qinit = qinit + make_qinit(
                x, y, Z,
                geojson,
                avid,
                AVAC["d0s"],
                AVAC.get("z0", 2000),
                AVAC.get("d0g", 5/100/100)
            )
    else:
        qinit = make_qinit(
            x, y, Z,
            geojson,
            avid,
            AVAC["d0s"],
            AVAC.get("z0", 2000),
            AVAC.get("d0g", 5/100/100)
        )
    print(f"Saving {projdir/AVAC['qinit']}", end="... ", flush=True)
    X, Y = np.meshgrid(x, y, copy=False)
    np.savetxt(projdir/AVAC['qinit'], np.column_stack((X.flatten(), Y.flatten(), qinit.flatten())))
    print("Saved.", flush=True)

    if plot:
        ext = x.min(), x.max(), y.min(), y.max()
        plt.imshow(Z, extent=ext)
        im = plt.imshow(np.where(qinit==0, np.nan, qinit), extent=ext, cmap=plt.cm.Reds)
        plt.xlabel("$x$ [m]")
        plt.ylabel("$y$ [m]")
        plt.title("Avalanche panels and depth")
        c = plt.colorbar(im)
        c.set_label("Snow fall over 3 days $d_0$ ($T=300$ y)")
        plt.show()


def make_qinit(x, y, Z, geojson, avid, d0s, z0, d0g):
    """From a table of polygons, describe qinit (depth, speeds are 0).

    Parameters
    ----------
    x: np.ndarray
        The x-coordinates of the bathymetry (1D)
    y: np.ndarray
        The y-coordinates of the bathymetry (1D)
    Z: np.ndarray
        The z-coordinates of the bathymetry (2D)
    geojson: np.ndarray
        A table of the polygons (avalanche_id, x, y)
    avid:
        The avalanche id of intrest (only one at a time)

    Returns
    ------
    d0: np.ndarray
        The 2D array of initial depth

    Reference
    ---------
    The empirical method used to compute the initial depth is described here:
    http://www.toraval.ch/articles/trad2.pdf
    """
    d0 = np.zeros((y.size, x.size), dtype=np.float32)
    ix, x_all, y_all = geojson
    ix = ix.astype(np.uint8)

    if avid not in ix:
        raise ValueError(f"Avalanche #{avid} is not in {np.unique(ix)}")
    # select the <avid> avalanche polygon
    xi = x_all[avid==ix]
    yi = y_all[avid==ix]
    # Find the bounding box of the contour (trying to save some resources)
    i0 = (xi.min() <= x).argmax()
    j0 = (y <= yi.max()).argmax()
    i1 = x.size-1 - (x <= xi.max())[::-1].argmax()
    j1 = y.size-1 - (yi.min() <= y)[::-1].argmax()
    X, Y = np.meshgrid(x[i0:i1+1], y[j0:j1+1])
    # Check wich points are inside the polygon
    inside = mPath(np.column_stack((xi, yi))).contains_points(
        np.column_stack((X.flatten(), Y.flatten()))
    ).reshape(X.shape)
    inside = isotropic_dilation(inside, 1)  # To avoid underestimating the volume
    # Compute d0
    phi = dip(
        x[i0:i1+1],
        y[j0:j1+1],
        Z[j0:j1+1, i0:i1+1]
    )[inside[1:-1, 1:-1]].mean()
    d = d0s - d0g/100/100 * (Z[j0:j1+1, i0:i1+1][inside] - z0)
    d0[j0:j1+1, i0:i1+1][inside] = (d*0.291/((np.sin(phi)-0.202*np.cos(phi))*np.cos(phi))).mean()
    # Compute avalanche volume to log it
    dx = x[2:] - x[:-2]
    dy = y[:-2] - y[2:]
    V = ((d0[1:-1, 1:-1]*dx).T*dy).sum()
    with open(projdir/"avac"/"qinit.log", "a") as file:
        print(f"{avid=}: {V=:.2e} {d0.mean() = :.2e} {np.rad2deg(phi)=:.2f}", flush=True, file=file)
    return d0


def dip(x, y, Z):
    """
    Compute the highest dip given a regular grid of the bathymetry.

    Parameters
    ----------
    x: np.ndarray
        x-coordinates (1D)
    y: np.ndarray
        y-coordinates (1D)
    Z: np.ndarray
        z-coordinates (2D)

    Returns
    -------
    dip: np.ndarray
        The highest dip at every point (2D)
    """
    dx = x[2:] - x[:-2]
    dy = y[:-2] - y[2:]
    D = (((Z[1:-1, 2:]-Z[1:-1, :-2])/dx)**2 + ((Z[2:, 1:-1]-Z[:-2, 1:-1]).T/dy).T**2)
    m = ~np.isclose(D, 0)
    D[m] = np.arccos(np.sqrt(D[m]) / np.sqrt(D[m]*(1+D[m])))
    D[~m] = 0
    return D


def _talweg(Z, visited, wmax: int, w=1, Zend=-float("inf"), tol_kw=dict()):
    """ Recursive utility to find the talweg. """
    x, y = visited[-1]
    x0 = max(0, x-w)
    y0 = max(0, y-w)
    x2 = min(Z.shape[1], x+w)
    y2 = min(Z.shape[0], y+w)
    Zs = Z[y0:y2+1, x0:x2+1]
    nx = np.argmin(Zs)
    ym = y0 + nx // Zs.shape[0]
    xm = x0 + nx % Zs.shape[0]
    if (xm, ym) in visited:
        w = w + 1
    else:
        w = 1
    if w >= wmax or (Z[ym, xm]<=Zend):
        return visited
    visited.append((xm, ym))
    return _talweg(Z, visited, wmax, w, Zend, tol_kw)


def talweg(Z, i, j, wstart=1, wmax=None, Zend=-float("inf"), rec_limit=int(1e5), tol_kw=None):
    """
    Find the talweg on bathymetry Z starting from the point (i,j).

    Parameters
    ----------
    Z: np.ndarray (2D)
    i: int
    j: int

    Returns
    -------
    visited: List[int, int]
        The coordinates of the points of the talweg
    """
    prev_rec_limit = getrecursionlimit()
    setrecursionlimit(rec_limit)
    if wmax is None:
        wmax = min(Z.shape)
    if tol_kw is None:
        tol_kw = dict(atol=0, rtol=0)
    visited = _talweg(Z, [(i, j)], Zend=Zend, w=wstart, wmax=wmax, tol_kw=tol_kw)
    setrecursionlimit(prev_rec_limit)
    return visited


def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("avid", nargs="?", default=1, type=int)
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    write_qinit(**args.__dict__)


if __name__ == "__main__":
    main()

