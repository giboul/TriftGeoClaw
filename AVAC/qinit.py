#!/usr/bin/env python
from pathlib import Path
from argparse import ArgumentParser
from yaml import safe_load
import numpy as np
from matplotlib.path import Path as mPath
from matplotlib import pyplot as plt
from skimage.morphology import isotropic_dilation


projdir = Path().absolute().parent
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    topoconfig = config["TOPM"]
    config = config["AVAC"]

def write_qinit(avid, plot=False):

    path = projdir / topoconfig["bathymetry"]
    print(f"\tINFO: Opening {path}... ", end="")
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
    z = np.loadtxt(path, skiprows=6).reshape(ny, nx)

    qinit = make_qinit(X, Y, indices=avid, plot=plot)
    print("Saving qinit.xyz...", end=" ", flush=True)
    np.savetxt("qinit.xyz", np.column_stack((X.flatten(), Y.flatten(), qinit.flatten())))
    print("Saved.")
    qinit[qinit <= 0] = float("nan")

    if plot:
        ext = xmin, x.max(), ymin, y.max()
        plt.imshow(z, extent=ext)
        plt.imshow(qinit, extent=ext, cmap=plt.cm.Blues)
        bounds = config.get("bounds") or topoconfig["bounds"]
        plt.scatter((bounds["xmin"], bounds["xmax"]), (bounds["ymin"], bounds["ymax"]))
        plt.legend()
        plt.show()


def make_qinit(X, Y, indices="", plot=False):
    Z = np.zeros_like(X, dtype=np.float16)
    print("Loading avalances.csv...", end=" ")
    ix, x_all, y_all = np.loadtxt(projdir / config["avalanches"])
    print("Loaded.")
    ix = ix.astype(np.uint8)
    if not indices:
        indices = np.unique(ix)
    else:
        indices = [int(indices)]
    for _i, i in enumerate(indices):
        print(f"Setting avalanche {i} ({_i+1}/{len(indices)})", end="\r")
        if i not in ix:
            raise ValueError(f"Avalanche #{i} is out o bounds {ix.min(), ix.max()}")
        x = x_all[i==ix]
        y = y_all[i==ix]
        path = mPath(np.column_stack((x, y)))
        inside = path.contains_points(np.column_stack((X.flatten(), Y.flatten())))
        inside = inside.reshape(X.shape)
        inside = isotropic_dilation(inside, 2)
        Z[inside] = 3
        if plot:
            plt.plot(x, y, label=i)
    print()
    return Z


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("avid", nargs="?", default="")
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    write_qinit(**args.__dict__)

