#!/usr/bin/env python
import json
from argparse import ArgumentParser
import AddSetrun
import numpy as np
from matplotlib.path import Path as mPath
from matplotlib import pyplot as plt
from skimage.morphology import isotropic_dilation
params = AddSetrun

def write_qinit(avid, plot=False):

    # Temporary directory
    ifile = "bathy_with_dam.asc"

    print(f"\tINFO: Opening {ifile}... ", end="")
    def readline(f, t):
        s = f.readline()
        return t(s[:s.find(" ")])

    with open(ifile) as file:
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
    z = np.loadtxt(ifile, skiprows=6).reshape(ny, nx)

    qinit = make_qinit(X, Y, indices=avid, plot=plot)
    # qinit[(X-2.668e6)**2 + (Y-1.1695e6)**2 <= 5e2**2] = 10
    # qinit[(X-2.671e6)**2 + (Y-1.1730e6)**2 <= 5e2**2] = 3
    # qinit[(X-2.671e6)**2 + (Y-1.1695e6)**2 <= 5e2**2] = 3
    # qinit[(X-2.672e6)**2 + (Y-1.1710e6)**2 <= 5e2**2] = 3
    np.savetxt("qinit.xyz", np.vstack((X.flatten(), Y.flatten(), qinit.flatten())).T)
    qinit[qinit <= 0] = float("nan")

    if plot is False:
        return None
    ext = xmin, x.max(), ymin, y.max()
    plt.imshow(z, extent=ext)
    plt.imshow(qinit, extent=ext, cmap=plt.cm.Blues)
    plt.scatter((list(params.bounds.values())[:2]),
                (list(params.bounds.values())[2:]))
    plt.show()


def make_qinit(X, Y, indices="", plot=False):
    Z = np.zeros_like(X, dtype=np.float16)
    print("Loading avalances.csv...", end=" ")
    ix, x_all, y_all = read_avalanches()
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
        path = mPath(np.vstack((x, y)).T)
        inside = path.contains_points(np.vstack((X.flatten(), Y.flatten())).T)
        inside = inside.reshape(X.shape)
        inside = isotropic_dilation(inside, 2)
        Z[inside] = 3
        if plot:
            plt.plot(x, y)
    print()
    return Z


def read_avalanches():
    with open("avalanches.geojson", "r") as file:
        data = json.load(file)
    
    coords = []
    for e in data["features"]:
        ix = e['properties']['id']
        points = e['geometry']['coordinates'][0][0]
        for x, y in points:
            coords.append([ix-1, x, y])
    
    return np.array(coords).T

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("avid", nargs="?", default="")
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    write_qinit(args.avid, args.plot)

