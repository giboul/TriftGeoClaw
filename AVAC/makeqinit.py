#!/usr/bin/env python
from argparse import ArgumentParser
import AddSetrun
import numpy as np
from matplotlib.path import Path
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

    qinit = np.zeros_like(X, dtype=np.float16)
    insert_avalanches(X, Y, qinit, indices=avid)
    # qinit[(X-2.668e6)**2 + (Y-1.1695e6)**2 <= 5e2**2] = 10
    # qinit[(X-2.671e6)**2 + (Y-1.1730e6)**2 <= 5e2**2] = 3
    # qinit[(X-2.671e6)**2 + (Y-1.1695e6)**2 <= 5e2**2] = 3
    # qinit[(X-2.672e6)**2 + (Y-1.1710e6)**2 <= 5e2**2] = 3
    np.savetxt("qinit.xyz", np.vstack((X.flatten(), Y.flatten(), qinit.flatten())).T)
    qinit[qinit <= 0] = float("nan")

    if plot == False:
        return None
    plt.figure(layout="tight")
    plt.imshow(z, extent=(xmin, x.max(), ymin, y.max()))
    plt.imshow(qinit, extent=(X.min(), X.max(), Y.min(), Y.max()))
    plt.scatter((list(params.bounds.values())[:2]),
                (list(params.bounds.values())[2:]))
    plt.show()


def insert_avalanches(X, Y, Z, indices=""):
    ix, x_all, y_all = np.loadtxt("avalanches.csv").T
    ix = ix.astype(np.uint8)
    if not indices:
        indices = np.unique(ix)
    else:
        indices = [int(indices)]
    for i in indices:
        if i not in ix:
            raise ValueError(f"Avalanche #{i} is out o bounds {ix.min(), ix.max()}")
        x = x_all[i==ix]
        y = y_all[i==ix]
        path = Path(np.vstack((x, y)).T)
        inside = path.contains_points(np.vstack((X.flatten(), Y.flatten())).T)
        inside = inside.reshape(X.shape)
        inside = isotropic_dilation(inside, 2)
        Z[inside] = 3
        # plt.imshow(inside, extent=(X.min(), X.max(), Y.min(), Y.max()))
        # plt.plot(x, y)
        # plt.show()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("avid", nargs="?", default="")
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    write_qinit(args.avid, args.plot)

