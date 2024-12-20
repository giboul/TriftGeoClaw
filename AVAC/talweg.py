from typing import List, Tuple, Iterable, Dict
from pathlib import Path
from sys import getrecursionlimit, setrecursionlimit
import numpy as np
from matplotlib import pyplot as plt

def _talweg(Z: Iterable[float], visited: List[Tuple], wmax: int,
            w: int=1, Zend: float=-float("inf"), tol_kw: Dict=dict()) -> List[Tuple]:
    x, y = visited[-1]
    x0 = max(0, x-w)
    y0 = max(0, y-w)
    x2 = min(Z.shape[1], x+w)
    y2 = min(Z.shape[0], y+w)
    Zs = Z[y0:y2+1, x0:x2+1]
    nx = np.argmin(Zs)
    ym = y0 + nx // Zs.shape[0]
    xm = x0 + nx % Zs.shape[0]
    # print(Z.shape, (f"Z[{y},{x}]={Z[y,x]}", f"Z[{ym}, {xm}]={Z[ym,xm]}"))
    if (xm, ym) in visited:
        w = w + 1
    else:
        w = 1
    if w >= wmax or (Z[ym, xm]<=Zend):
        return visited
    visited.append((xm, ym))
    return _talweg(Z, visited, wmax, w, Zend, tol_kw)

def talweg(Z: Iterable[float], i: int, j: int,
           wstart: int=1, wmax: int=None,
           rec_limit:int=int(1e5), tol_kw=None):
    prev_rec_limit = getrecursionlimit()
    setrecursionlimit(rec_limit)
    if wmax is None:
        wmax = min(Z.shape)
    if tol_kw is None:
        tol_kw = dict(atol=0, rtol=0)
    visited = _talweg(Z, [(i, j)], Zend=1775, w=wstart, wmax=wmax, tol_kw=tol_kw)
    setrecursionlimit(prev_rec_limit)
    return visited

path = Path(__file__).parent.parent / "TOPM" / "bathymetry.asc"
print(f"INFO: Opening {path}... ", end="")
def readline(f, t):
    s = f.readline()
    return t(s[:s.find(" ")])
with open(path) as file:
    nx = readline(file, int)
    ny = readline(file, int)
    xmin = readline(file, float)
    ymin = readline(file, float)
    resolution = readline(file, float)
    nodatavalue = readline(file, float)
print("Loaded.")
Z = np.loadtxt(path, skiprows=6).reshape(ny, nx)

pts = talweg(Z, 6*Z.shape[0]//10, 6*Z.shape[1]//10, wmax=131)
# print(pts)

plt.imshow(Z)
plt.plot(*list(zip(*pts)), c='r')
plt.show()

