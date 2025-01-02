from typing import List, Tuple, Iterable, Dict
from pathlib import Path
from sys import getrecursionlimit, setrecursionlimit
import numpy as np
from matplotlib import pyplot as plt
plt.style.use("seaborn-v0_8-paper")

def _talweg(Z: Iterable[float], visited: List[Tuple], wstart: int=1, wmax: int=np.inf,
            w: int=1, Zend: float=-float("inf"), tol_kw: Dict=dict()) -> List[Tuple]:
    x, y = visited[-1]
    x0 = max(0, x-w)
    y0 = max(0, y-w)
    x2 = min(Z.shape[1]-1, x+w)
    y2 = min(Z.shape[0]-1, y+w)
    Zs = Z[y0:y2+1, x0:x2+1]
    nx = np.argmin(Zs)
    ym = y0 + nx // Zs.shape[0]
    xm = x0 + nx % Zs.shape[0]
    # print(Z.shape, (f"Z[{y},{x}]={Z[y,x]}", f"Z[{ym}, {xm}]={Z[ym,xm]}"))
    if (xm, ym) in visited:
        w = w + 1
    else:
        visited.append((xm, ym))
        w = wstart
    if w >= wmax or (Z[ym, xm]<=Zend):
        return visited
    return _talweg(Z, visited, wstart, wmax, w, Zend, tol_kw)

def talweg(Z: Iterable[float], i: int, j: int,
           wstart: int=1, wmax: int=None,
           rec_limit:int=int(1e5), tol_kw=None, Zend=-np.inf):
    prev_rec_limit = getrecursionlimit()
    setrecursionlimit(rec_limit)
    if wmax is None:
        wmax = min(Z.shape)
    if tol_kw is None:
        tol_kw = dict(atol=0, rtol=0)
    visited = _talweg(Z, [(i, j)], Zend=Zend, wstart=wstart, w=wstart, wmax=wmax, tol_kw=tol_kw)
    setrecursionlimit(prev_rec_limit)
    return np.array(visited).T

def main():
    projdir = Path(__file__).parent.parent
    path =  projdir / "natural_bathymetry.bin"
    print(f"INFO: Opening {path}... ", end="")
    def readline(f, t):
        s = f.readline()
        return t(s[:s.find(" ")])
    with open(path.with_suffix(".data")) as file:
        nx = readline(file, int)
        ny = readline(file, int)
        xmin = readline(file, float)
        ymin = readline(file, float)
        resolution = readline(file, float)
        nodatavalue = readline(file, float)
    print("Loaded.")
    ext = xmin, xmin+nx*resolution, ymin, ymin+ny*resolution
    Z = np.fromfile(path, dtype=np.float64).reshape(ny, nx, order="F")
    print(f"{Z.max() = }")
    print(f"{Z.min() = }")

    geojson = np.loadtxt(projdir / "TOPM" / "avalanches.csv")

    fig, ax = plt.subplots(nrows=3, layout="tight", figsize=(6, 11))
    for ix in np.unique(geojson[0]):
        xa, ya = geojson[1:, geojson[0, :] == ix]
        xi, yi = talweg(Z,
                        int((xa.mean()-xmin)/resolution),
                        int(ny-1-(ya.mean()-ymin)/resolution),
                        wstart=10,
                        Zend=1767)
        x = xmin + xi*resolution
        y = ymin + (ny-1-yi)*resolution
        # print(pts)
        ax[0].fill(*geojson[:, geojson[0]==ix][1:])
        ax[0].imshow(Z, extent=ext)
        ax[0].plot(x, y)
        dl = np.sqrt(np.diff(x)**2+np.diff(y)**2)
        l = np.hstack((0, np.cumsum(dl)))
        z = Z[yi, xi]
        ax[1].plot(l, z, '-', ms=5, label=int(ix))
        ax[2].plot((l[1:]+l[:-1])/2/l[-1], -np.rad2deg(np.arctan((z[1:]-z[:-1])/dl)))
        # ax[2].plot((l[1:]+l[:-1])/2, (z[1:]-z[:-1])/dl)
    ax[1].set_aspect("equal")
    ax[1].legend(ncols=5, title="Identifiant de l'avalanche")
    ax[1].set_xlabel(r"Distance parcourue $\ell$ [m]")
    ax[1].set_ylabel("Altitude $z$ [m.s.m.]")
    ax[1].set_ylim(1767, None)
    ax[1].margins(x=0, y=0)
    ax[2].set_xlabel(r"Distance paurcourue normalisée $\ell/\ell_f$ [m]")
    ax[2].set_ylabel(r"Pente $\psi$ [°]")
    ax[2].axline((0, 27.5), slope=0, ls='-.', c='grey')
    plt.show()

if __name__ == "__main__":
    main()
