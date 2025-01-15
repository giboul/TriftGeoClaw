#!/usr/bin/env python
from sys import getrecursionlimit, setrecursionlimit
from pathlib import Path
from argparse import ArgumentParser
from yaml import safe_load
import numpy as np
from matplotlib.path import Path as mPath
from matplotlib import pyplot as plt
from skimage.morphology import isotropic_dilation


projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    TOPM = config["TOPM"]
    AVAC = config["AVAC"]

def write_qinit(avid="", plot=False):

    path = projdir / "TOPM" / "natural_bathymetry.bin"
    print(f"\tINFO: Opening {path}... ", end="", flush=True)
    def readline(f, t):
        s = f.readline()
        ix = s.find(" ")
        return t(s[:ix]), s[ix:]
    with open(path.with_suffix(".data")) as file:
        nx, _ = readline(file, int)
        ny, _ = readline(file, int)
        xmin, _ = readline(file, float)
        ymin, _ = readline(file, float)
        resolution, _ = readline(file, float)
        nodatavalue, _ = readline(file, float)
    x = xmin + np.arange(nx)*resolution
    y = ymin + np.arange(ny)[::-1]*resolution
    Z = np.fromfile(path, dtype=np.float16).reshape(ny, nx, order="F").astype(np.float32)
    print("Loaded.", flush=True)

    # for p in (projdir/"AVAC").glob("qinit*.xyz"):
    #     p.unlink()
    geojson = np.loadtxt(projdir / "TOPM" / "avalanches.csv")
    if avid:
        avids = [int(a) for a in avid.split(",")]
    else:
        avids = np.unique(geojson[0].astype(np.uint8))
    filename = projdir/"AVAC"/f"qinit{avid}.xyz"
    qinit = make_qinit(x, y, Z, geojson, indices=avids)
    if plot:
        ext = x.min(), x.max(), y.min(), y.max()
        plt.figure(layout="tight")
        plt.imshow(Z, extent=ext)
        # for ix in avids:
        #     x, y = geojson[:, geojson[0, :] == ix][1:, :].mean(axis=1)
        #     iy, ix = np.unravel_index(((X-x)**2+(Y-y)**2).argmin(), shape=Z.shape)
        #     _t0 = pc()
        #     xt, yt = np.array(talweg(Z, ix, iy, Zend=TOPM["lake_alt"])).T
        #     # print(f"{pc() - _t0}")
        #     xt = xmin + xt*resolution
        #     yt = ymin+ny*resolution - yt*resolution
        #     l, = plt.plot(xt, yt, c='r')
        # im = plt.imshow(np.where(qinit==0, np.nan, qinit), extent=ext, cmap=plt.cm.Reds, zorder=l.get_zorder()+1)
        im = plt.imshow(np.where(qinit==0, np.nan, qinit), extent=ext, cmap=plt.cm.Reds)
        plt.xlabel("$x$ [m]")
        plt.ylabel("$y$ [m]")
        plt.title("Avalanche panels and depth")
        c = plt.colorbar(im)
        c.set_label("Snow fall over 3 days $d_0$ ($T=300y$)")
        # plt.legend()
        plt.show()
    print(f"\tINFO: Saving {filename}", end="... ", flush=True)
    X, Y = np.meshgrid(x, y, copy=False)
    np.savetxt(filename, np.column_stack((X.flatten(), Y.flatten(), qinit.flatten())), fmt="%.9e")
    print("Saved.", flush=True)


def make_qinit(x, y, Z, geojson, indices):

    dx = (x[2:] - x[:-2])[0]
    dy = (y[:-2] - y[2:])[0]
    H = np.zeros((y.size, x.size), dtype=np.float32)
    print("\tINFO: Loading avalances.csv...", end=" ", flush=True)
    ix, x_all, y_all = geojson
    print("Loaded.", flush=True)
    ix = ix.astype(np.uint8)

    X, Y = np.meshgrid(x, y)
    for _i, i in enumerate(indices):
        # print(f"\tINFO: Setting avalanche {i} ({_i+1}/{len(indices)})", end="...\r", flush=True)
        if i not in ix:
            raise ValueError(f"Avalanche #{i} is not in {np.unique(ix)}")
        xi = x_all[i==ix]
        yi = y_all[i==ix]
        i0 = (xi.min() <= x).argmax()
        j0 = (y <= yi.max()).argmax()
        i1 = x.size-1 - (x <= xi.max())[::-1].argmax()
        j1 = y.size-1 - (yi.min() <= y)[::-1].argmax()
        X, Y = np.meshgrid(x[i0:i1+1], y[j0:j1+1])
        inside = mPath(np.column_stack((xi, yi))).contains_points(
            np.column_stack((X.flatten(), Y.flatten()))
        ).reshape(X.shape)
        inside = isotropic_dilation(inside, 2)
        p = dip(dx, dy, Z[j0:j1+1, i0:i1+1])[inside[1:-1, 1:-1]].mean()
        d = 2.0 - 5/100/100 * (Z[j0:j1+1, i0:i1+1][inside] - 2000)  # T = 300 years, Western Bernese Oberland
        H[j0:j1+1, i0:i1+1][inside] = (d*0.291/((np.sin(p)-0.202*np.cos(p))*np.cos(p))).mean()
        if i == 25 or i == 24:  # ice avalanches
            H[j0:j1+1, i0:i1+1][inside] = 15 * 900./300.
        d0 = H[j0:j1+1, i0:i1+1][inside]
        V = (d0*dx*dy).sum()
        with open(projdir/"AVAC"/"qinit.log", "a") as file:
            print(f"\t\t{i=}: {V=:.2e} {d0.mean() = :.2e} {np.rad2deg(p)=:.2f}", flush=True, file=file)
    print()
    return H

def dip(dx, dy, Z):
    Z = (((Z[1:-1, 2:]-Z[1:-1, :-2])/dx)**2 + ((Z[2:, 1:-1]-Z[:-2, 1:-1])/dy)**2)
    m = ~np.isclose(Z, 0)
    Z[m] = np.arccos(np.sqrt(Z[m]) / np.sqrt(Z[m]*(1+Z[m])))
    Z[~m] = 0
    return Z

def _talweg(Z, visited, wmax: int, w=1, Zend=-float("inf"), tol_kw=dict()):
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
    prev_rec_limit = getrecursionlimit()
    setrecursionlimit(rec_limit)
    if wmax is None:
        wmax = min(Z.shape)
    if tol_kw is None:
        tol_kw = dict(atol=0, rtol=0)
    visited = _talweg(Z, [(i, j)], Zend=Zend, w=wstart, wmax=wmax, tol_kw=tol_kw)
    setrecursionlimit(prev_rec_limit)
    return visited


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("avid", nargs="?", default="")
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    write_qinit(**args.__dict__)

