import numpy as np
from matplotlib import pyplot as plt


Lf = 3
sg = 0.4
ls = 0.85
l0 = 0.3
h0 = 0.2
theta = np.pi/6
slope = np.tan(theta)
H = np.sin(theta) * ls + h0
L = H / np.tan(theta)
lp = np.sqrt(H**2 + L**2)

step = 0.05
X, Y = np.meshgrid(np.arange(0, Lf+step, step=step), np.arange(0, 0.2+step, step=step)[::-1])

Z = np.zeros_like(X, dtype=np.float64)
Z[X <= L] = H - X[X <= L] * slope

H = np.maximum(Z, h0) - Z
H = np.zeros_like(Z)
H[X <= l0] = sg * (1-l0*slope/2/sg) / np.cos(theta)


def main():
    fig, ax = plt.subplots(ncols=2)
    ax[0].set_aspect("equal")
    ax[0].plot(X[0, :], Z[0, :])
    ax[0].plot(X[0, :], Z[0, :] + H[0, :])
    ax[1].imshow(Z)
    plt.show()

    x, y, z, h = [a.flatten() for a in (X, Y, Z, H)]
    asc_header = "\n".join((
        f"{X.shape[1]} ncols",
        f"{Y.shape[0]} nrows",
        f"{x.min()} xllcenter",
        f"{y.min()} yllcenter",
        f"{step} cellsize",
        f"{999999} nodata_value"
    ))

    np.savetxt("bathymetry.asc", z, header=asc_header, comments="")
    np.savetxt("qinit.xyz", np.column_stack((x, y, h)))


if __name__ == "__main__":
    main()

