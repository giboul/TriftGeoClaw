from pathlib import Path
from argparse import ArgumentParser
from yaml import safe_load
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.path import Path as mPath
from clawpack.visclaw import gridtools
from clawpack.pyclaw import solution


parser = ArgumentParser()
parser.add_argument("avid", type=str, nargs="?", default="")
avid = parser.parse_args().avid

projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    AVAC = config["AVAC"]
    TOPM = config["TOPM"]
    TSUL = config["TSUL"]

contour = np.loadtxt(projdir / "TOPM" / "contour2.xy").T
contour_lake = np.loadtxt(projdir / "TOPM" / "contour1.xy").T

outdir = f"_output{avid}"

files = sorted(list(Path(projdir/"AVAC"/outdir).glob("fort.t*")))
timesAVAC = np.zeros(len(files), np.float64)
for i, file in enumerate(files):
    with open(file) as file:
        t = float([e for e in file.readlines()[0].split() if e != ""][0])
        timesAVAC[i] = t

times = np.loadtxt(projdir/"TSUL"/outdir/"times.txt")

with open(projdir/"TSUL"/"claw.data") as file:
    lines = file.readlines()[5:]
clawdata_trans = dict(T=True, F=False)
clawdata = []
for line in lines:
    elements = line[:(line.find("=")+1 or 0)-1].split()
    if elements:
        clawdata.append([])
        for e in elements:
            out = clawdata_trans.get(e)
            if out is None:
                out = eval(e)
            clawdata[-1].append(out)
ndim, (xmin, ymin), (xmax, ymax), (numx, numy), *_ = clawdata

x = np.linspace(xmin, xmax, 100)
y = np.linspace(ymin, ymax, 100)
X, Y = np.meshgrid(x, y)
dx = (xmax-xmin)/numx
dy = (ymax-ymin)/numy
dx = (x[2:] - x[:-2])/2
dy = (y[2:] - y[:-2])/2

def read_TSUL(i, outdir=projdir/"TSUL"/outdir):

    frame = solution.Solution(i, path=outdir, file_format=TSUL["out_format"])
    q = gridtools.grid_output_2d(frame, lambda x: x, X, Y, levels='all',return_ma=True)

    return q


def interp_AVAC_contour(t, x, y, outdir=projdir/"AVAC"/outdir):
    it = np.clip((timesAVAC <= t).sum()-1, 0, timesAVAC.size-2)
    frame_sol = solution.Solution(it, path=outdir, file_format=AVAC["out_format"])
    q1 = gridtools.grid_output_2d(
        frame_sol,
        lambda q: q,
        x, y,
        levels = "all",
        return_ma=True
    )
    frame_sol = solution.Solution(it+1, path=outdir, file_format=AVAC["out_format"])
    q2 = gridtools.grid_output_2d(
        frame_sol,
        lambda q: q,
        x, y,
        levels = "all",
        return_ma=True
    )
    t = timesAVAC[it] + (t - timesAVAC[it])/(timesAVAC[it+1] - timesAVAC[it])
    q2 = q1 + (q2 - q1) * (t - timesAVAC[it])/(timesAVAC[it+1] - timesAVAC[it])
    # fig, axes = plt.subplots(ncols=q2.shape[0])
    # for ax, var in zip(axes, q2):
    #     ax.plot(var)
    # fig.show()

    return q2


def read_TSUL_contour(i, x, y, outdir=projdir/"TSUL"/outdir):
    frame_sol = solution.Solution(i, path=outdir, file_format=TSUL["out_format"])
    q = gridtools.grid_output_2d(
        frame_sol,
        lambda q: q,
        x, y,
        levels = "all",
        return_ma=True
    )
    return q


def divide(a, b, fill=0.):
    c = np.empty(a.shape)
    m = np.isclose(b, 0.)
    c[m] = fill
    m = ~m
    c[m] = a[m] / b[m]
    return c


nsols = len(list(Path(projdir/"TSUL"/outdir).glob("fort.t*")))
solutions = [read_TSUL(i) for i in range(nsols)]

g = 9.81
rhow = 1000
rhos = 300
Elac = np.zeros(times.size, np.float64)
Etav = np.zeros(times.size, np.float64)
Eavc = np.zeros(times.size, np.float64)
lake = (solutions[0][3] <= TOPM["lake_alt"]) | np.isclose(solutions[0][3], TOPM["lake_alt"])
lake = mPath(contour_lake.T).contains_points(np.column_stack((X.flatten(), Y.flatten()))).reshape(X.shape)
# lake = np.full_like(lake, True)
xc, yc = contour
dxc  = np.hstack((xc[1]-xc[-1], xc[2:] - xc[:-2], xc[0]-xc[-2]))
dyc  = np.hstack((yc[1]-yc[-1], yc[2:] - yc[:-2], yc[0]-yc[-2]))
dl = np.sqrt(dxc**2 + dyc**2)
nx, ny = divide(-dyc, dl), divide(dxc, dl)

def trapint(x, y):
    ym = (y[1:] + y[:-1])/2
    xm = (x[1:] + x[:-1])/2
    return xm, np.cumsum(ym/2 * np.diff(x))

def intdS(f):
    f = f[1:-1, 1:-1]
    f = (f*dx).T
    f = (f*dy).T
    return f[lake[1:-1, 1:-1]].sum()

def intdl(f):
    return (f*dl).sum()

for i, (q, t) in enumerate(zip(solutions, times)):
    fig, (ax1, ax2) = plt.subplots(ncols=2)
    print(i, t)
    # Lake
    h, hu, hv, s = q
    eta = s - TOPM["lake_alt"]
    eta = np.where(~lake, np.nan, eta)
    u = divide(hu, h)
    v = divide(hv, h)
    u2 = u**2+v**2
    # print(np.nanmin(eta), np.nanmax(eta))
    ax1.imshow(eta, origin="lower", extent=(X.min(), X.max(), Y.min(), Y.max()))
    Elac[i] = 1/2*rhow*intdS(g*eta**2 + h*u2)
    # Avalanche TSUL
    h, hu, hv, s = read_TSUL_contour(i, xc, yc)
    eta = s - TOPM["lake_alt"]
    u = divide(hu, h)
    v = divide(hv, h)
    un = u*nx+v*ny
    u2 = u**2+v**2
    ax2.plot(un, ls="-", marker='o', ms=2, mfc="w", label="TSUL")
    Etav[i] = 1/2*rhow*intdl(h*g*eta*un + h*u2*un)
    # Avalanche AVAC
    h, hu, hv, s = interp_AVAC_contour(t, xc, yc)
    eta = s - TOPM["lake_alt"]
    u = divide(hu, h)
    v = divide(hv, h)
    un = u*nx+v*ny
    u2 = u**2+v**2
    ax2.plot(un, ls="-", marker='o', ms=2, mfc="w", label="AVAC")
    Eavc[i] = 1/2*rhos*intdl(h*g*eta*un + h*u2*un)
    ax2.legend()
    plt.colorbar(ax1.scatter(xc, yc, c=u2))
    plt.show(block=False)
    fig.show()

tm, Eavc = trapint(times, Eavc)
tm, Etav = trapint(times, Etav)

fig, ax = plt.subplots()
ax.plot(tm, Eavc, '-o', label=r"$E_\text{avalanche}^\text{AVAC}$")
ax.plot(tm, Etav, '-.', label=r"$E_\text{avalanche}^\text{TSUL}$")
ax.plot(times, Elac, '-o', label=r"$E_\text{lake}$")
ax.legend()
plt.show()
