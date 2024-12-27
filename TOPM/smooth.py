import numpy as np


def cartesian_grid_interp(xll, yll, nx, ny, cs, Zt, x, y):
    xur = xll + nx*cs
    yur = yll + ny*cs
    w = np.clip(0, nx-2, ((x-xll)/(xur-xll)*(nx-1)).astype(np.int64))
    s = np.clip(0, ny-2, ((y-yll)/(yur-yll)*(ny-1)).astype(np.int64))
    xw = xll + (xur-xll)*w/(nx-1)
    ys = yll + (yur-yll)*s/(ny-1)
    xe = xw + (xur-xll)/(nx-1)
    yn = ys + (yur-yll)/(ny-1)
    Z = (
        + (x-xw)*((y-ys)*Zt[:, w+1][s+1, :].T).T
        + (x-xw)*((yn-y)*Zt[:, w+1][s  , :].T).T
        + (xe-x)*((y-ys)*Zt[:, w  ][s+1, :].T).T
        + (xe-x)*((yn-y)*Zt[:, w  ][s  , :].T).T
    ) / cs**2
    return Z

def conv2d(a, f):
    s = f.shape + tuple(np.subtract(a.shape, f.shape) + 1)
    strd = np.lib.stride_tricks.as_strided
    subM = strd(a, shape = s, strides = a.strides * 2)
    return np.einsum('ij,ijkl->kl', f, subM)


def smooth_cartesian_grid(xll, yll, nx, ny, cs, Zt, new_res):
    x = np.arange(xll, xll+(nx+1)*cs, step=new_res)
    y = np.arange(yll, yll+(ny+1)*cs, step=new_res)
    Z = cartesian_grid_interp(xll, yll, nx, ny, cs, Zt, x, y)
    x = xll + np.arange(nx+1)*cs
    y = yll + np.arange(ny+1)*cs
    Z = cartesian_grid_interp(xll, yll, nx, ny, cs, Zt, x, y)
    return Z


def main():
    from pathlib import Path
    from argparse import ArgumentParser

    topo = Path(__file__).parent / "bathymetry.asc"
    def readline(file):
        line = file.readline()
        line = line[:len(line)-1-line[::-1].find(" ")]
        return line.strip(" ")
    with open(topo) as file:
        nx = int(readline(file))
        ny = int(readline(file))
        xll = float(readline(file))
        yll = float(readline(file))
        cs = float(readline(file))
        ndv = float(readline(file))

    xur = xll+nx*cs
    yur = yll+ny*cs
    extent = xll, xur, yll, yur
    Z = np.loadtxt(topo, skiprows=6).reshape(ny, nx)

    # xi = np.linspace(xll, xll+nx*cs, num=5000, endpoint=True)
    # yi = np.linspace(yll, yll+ny*cs, num=5000, endpoint=True)
    # Zi = cartesian_grid_interp(xll, yll, nx, ny, cs, Z, xi, yi)
    kernel = np.array((
        (4, 4, 4, 4),
        (4, 2, 2, 4),
        (4, 2, 2, 4),
        (4, 4, 4, 4)
    ))
    kernel = kernel / kernel.sum()
    Zc = conv2d(Z, kernel)
    print(Z.shape, Zc.shape)

    parser = ArgumentParser()
    parser.add_argument("--plot", "-p", action="store_true")
    args = parser.parse_args()
    if args.plot is True:
        from matplotlib import pyplot as plt
        fig, (ax1, ax2) = plt.subplots(ncols=2, sharex=True, sharey=True)
        ax1.set_title("Original input")
        ax1.imshow(Z, extent=extent, cmap=plt.cm.RdBu)
        ax2.set_title("Interpolation output")
        ax2.imshow(Zc, extent=extent, cmap=plt.cm.RdBu)
        plt.show()


if __name__ == "__main__":
    main()
