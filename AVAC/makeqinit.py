#!/usr/bin/env python
import AddSetrun
import numpy as np
from matplotlib.path import Path
from matplotlib import pyplot as plt
params = AddSetrun

def write_qinit():

    # Temporary directory
    ifile = "bathy_with_dam.asc"

    print(f"\tINFO: Opening {ifile}... ")
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
    x = xmin + np.arange(nx)*resolution
    y = ymin + np.arange(ny)[::-1]*resolution
    X, Y = np.meshgrid(x, y)
    z = np.loadtxt(ifile, skiprows=6).reshape(ny, nx)


    qinit = np.zeros_like(X, dtype=np.float16)
    insert_avalanches(X, Y, qinit)
    # qinit[(X-2.668e6)**2 + (Y-1.1695e6)**2 <= 5e2**2] = 10
    # qinit[(X-2.671e6)**2 + (Y-1.1730e6)**2 <= 5e2**2] = 3
    # qinit[(X-2.671e6)**2 + (Y-1.1695e6)**2 <= 5e2**2] = 3
    # qinit[(X-2.672e6)**2 + (Y-1.1710e6)**2 <= 5e2**2] = 3
    np.savetxt("qinit.xyz", np.vstack((X.flatten(), Y.flatten(), qinit.flatten())).T)
    qinit[qinit <= 0] = float("nan")

    plt.imshow(z, extent=(xmin, x.max(), ymin, y.max()))
    plt.imshow(qinit, extent=(X.min(), X.max(), Y.min(), Y.max()))
    plt.scatter((list(params.bounds.values())[:2]),
                (list(params.bounds.values())[2:]))
    plt.show()


def insert_avalanches(X, Y, Z, indices="all"):
    ix, x_all, y_all = np.loadtxt("avalanches.csv").T
    if indices == "all":
        indices = ix
    for i in indices:
        x = x_all[i==ix]
        y = y_all[i==ix]
        path = Path(np.vstack((x, y)).T)
        inside = path.contains_points(np.vstack((X.flatten(), Y.flatten())).T)
        inside = inside.reshape(X.shape)
        Z[inside] = 3


if __name__ == "__main__":
    write_qinit()

