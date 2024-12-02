#!/usr/bin/env python
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

    path = projdir / TOPM["bathymetry"]
    print(f"INFO: Opening {path}... ", end="")
    def readline(f, t):
        s = f.readline()
        return t(s[:s.find(" ")])
    with open(path) as file:
        nx = readline(file, int)
        ny = readline(file, int)
        xmin = readline(file, float)
        ymin = readline(file, float)
        resolution = readline(file, int)
        nodatavalue = readline(file, float)
    print("Loaded.")
    x = xmin + np.arange(nx)*resolution
    y = ymin + np.arange(ny)[::-1]*resolution
    X, Y = np.meshgrid(x, y)
    Z = np.loadtxt(path, skiprows=6).reshape(ny, nx)

    for p in (projdir/"AVAC").glob("qinit*.xyz"):
        p.unlink()
    geojson = np.loadtxt(projdir / AVAC["avalanches"])
    if avid:
        avids = [int(a) for a in avid.split(",")]
    else:
        avids = np.unique(geojson[0].astype(np.int64))
    filename = projdir/"AVAC"/f"qinit{avid}.xyz"
    qinit = make_qinit(X, Y, Z, geojson, indices=avids)
    np.savetxt(filename, np.column_stack((X.flatten(), Y.flatten(), qinit.flatten())))
    if plot:
        ext = X.min(), X.max(), Y.min(), Y.max()
        plt.figure(layout="tight")
        plt.imshow(Z, extent=ext)
        im = plt.imshow(np.ma.MaskedArray(qinit, mask=qinit==0), extent=ext, cmap=plt.cm.Reds)
        plt.xlabel("$x$ [m]")
        plt.ylabel("$y$ [m]")
        plt.title("Avalanche panels and depth")
        c = plt.colorbar(im)
        c.set_label("Snow fall over 3 days $d_0$ ($T=300y$)")
        # plt.legend()
        plt.show()


def make_qinit(X, Y, Z, geojson, indices):

    H = np.zeros_like(X, dtype=np.float16)
    print("Loading avalances.csv...", end=" ")
    ix, x_all, y_all = geojson
    print("Loaded.")
    ix = ix.astype(np.uint8)
    d0s = 2.0 - 5/100/100 * (Z[1:-1, 1:-1] - 2000)  # T = 300 years, Western Bernese Oberland
    psi = dip(X, Y, Z)
    # _fig, _ax = plt.subplots()
    # _ax.set_axis_off()
    # _im = _ax.imshow(psi)
    # _ax.set_axis_off()
    # # plt.colorbar(_im)
    # _fig.savefig("/home/axel/Downloads/dip.pdf", bbox_inches="tight")
    # plt.show()
    # exit()

    for _i, i in enumerate(indices):
        print(f"Setting avalanche {i} ({_i+1}/{len(indices)})", end="\r")
        if i not in ix:
            raise ValueError(f"Avalanche #{i} is not in {np.unique(ix)}")
        x = x_all[i==ix]
        y = y_all[i==ix]
        path = mPath(np.column_stack((x, y)))
        inside = path.contains_points(np.column_stack((X.flatten(), Y.flatten())))
        inside = inside.reshape(X.shape)
        inside = isotropic_dilation(inside, 2)
        # H[1:-1, 1:-1][inside] = np.maximum(Z[1:-1, 1:-1][inside] - 1000, 0)*30/100/100
        p = psi[inside[1:-1, 1:-1]].mean()
        d = d0s[inside[1:-1, 1:-1]].mean()
        H[inside] = d * 0.291/(np.sin(p)-0.202*np.cos(p))
    print()
    return H

def dip(X, Y, Z):
    dX = X[1:-1, 2:] - X[1:-1,:-2]
    dY = Y[2:, 1:-1] - Y[:-2, 1:-1]
    gradZ = np.zeros((3, Z.shape[0]-2, Z.shape[1]-2), dtype=np.float64)
    gradZ[0] = (Z[1:-1, 2:]-Z[1:-1, :-2])/dX
    gradZ[1] = (Z[2:, 1:-1]-Z[:-2, 1:-1])/dY
    gradZ[2] = (gradZ[:2]**2).sum(axis=0)
    dip = np.zeros((Z.shape[0]-2, Z.shape[1]-2), dtype=np.float64)
    m = ~np.isclose(gradZ[2], 0)
    # dip[m] = np.arccos(np.linalg.norm(gradZ[:2], axis=0)[m] / np.linalg.norm(gradZ, axis=0)[m])
    dip[m] = np.arccos(np.sqrt(gradZ[2])[m] / np.sqrt(gradZ[2]*(1+gradZ[2]))[m])
    return dip


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("avid", nargs="?", default="")
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    write_qinit(**args.__dict__)

