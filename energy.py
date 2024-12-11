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

projdir = Path(__file__).parent
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    AVAC = config["AVAC"]
    TOPM = config["TOPM"]
    TSUL = config["TSUL"]

contour1 = np.loadtxt(projdir / "TOPM" / "contour1.xy").T
contour2 = np.loadtxt(projdir / "TOPM" / "contour2.xy").T

outdir = f"_output{avid}"

files = sorted(list(Path(projdir/"AVAC"/outdir).glob("fort.t*")))
timesAVAC = np.zeros(len(files), np.float64)
for i, file in enumerate(files):
    with open(file) as file:
        t = float([e for e in file.readlines()[0].split() if e != ""][0])
        timesAVAC[i] = t

times = np.loadtxt(projdir/"TSUL"/outdir/"times.txt")

def read_clawdata(outdir=projdir/"TSUL"/f"_output{avid}"):
    clawdata_trans = dict(T=True, F=False)
    clawdata = dict()
    with open(outdir/"claw.data") as file:
        lines = [l for l in file.readlines() if "=" in l]
        for line in lines:
            value, key = line.split("=: ")
            key = key[:key.find(" ")]
            value = [v for v in value.split() if v]
            for e, element in enumerate(value):
                try:
                    value[e] = eval(element)
                except Exception:
                    value[e] = clawdata_trans[element]
            clawdata[key] = value[0] if len(value)==1 else value
    return clawdata
clawdata = read_clawdata()
for k, v in clawdata.items():
    print(f"{k}: {v}")
xmin, ymin = clawdata["lower"]
xmax, ymax = clawdata["upper"]
numx, numy = clawdata["num_cells"]

x = np.linspace(xmin, xmax, 100)
y = np.linspace(ymin, ymax, 100)
X, Y = np.meshgrid(x, y)
extent = X.min(), X.max(), Y.min(), Y.max()
dx = (xmax-xmin)/numx
dy = (ymax-ymin)/numy
dx = (x[2:] - x[:-2])/2
dy = (y[2:] - y[:-2])/2

def read_TSUL(i, outdir=projdir/"TSUL"/outdir):
    frame = solution.Solution(int(i), path=outdir, file_format=TSUL["out_format"])
    return gridtools.grid_output_2d(frame, lambda x: x, X, Y, levels='all',return_ma=True)


def interp_AVAC_contour(t, x, y, outdir=projdir/"AVAC"/outdir):
    it = np.clip((timesAVAC <= t).sum()-1, 0, timesAVAC.size-2)
    frame_sol = solution.Solution(int(it), path=outdir, file_format=AVAC["out_format"])
    q1 = gridtools.grid_output_2d(
        frame_sol,
        lambda q: q,
        x, y,
        levels = "all",
        return_ma=True
    )
    frame_sol = solution.Solution(int(it+1), path=outdir, file_format=AVAC["out_format"])
    q2 = gridtools.grid_output_2d(
        frame_sol,
        lambda q: q,
        x, y,
        levels = "all",
        return_ma=True
    )
    q2 = q1 + (q2 - q1) * (t - timesAVAC[it])/(timesAVAC[it+1] - timesAVAC[it])
    # fig, axes = plt.subplots(ncols=q2.shape[0])
    # for ax, var in zip(axes, q2):
    #     ax.plot(var)
    # fig.show()

    return q2


def read_TSUL_contour(i, x, y, outdir=projdir/"TSUL"/outdir):
    frame_sol = solution.Solution(int(i), path=outdir, file_format=TSUL["out_format"])
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
Elake = np.zeros(times.size, np.float64)
Etsul = np.zeros(times.size, np.float64)
Eavac = np.zeros(times.size, np.float64)
Vavc = np.zeros(times.size, np.float64)
Vtav = np.zeros(times.size, np.float64)
Vlak = np.zeros(times.size, np.float64)

h0, _, _, eta0 = solutions[0]
z = eta0 - h0
lake = h0 >= 1e-5
lake = mPath(contour1.T).contains_points(np.column_stack((X.flatten(), Y.flatten()))).reshape(X.shape)
xc, yc = contour2
inside = (extent[0] <= xc) & (xc <= extent[1]) & (extent[2] <= yc) & (yc <= extent[3])
xc = xc[inside]
yc = yc[inside]
dxc  = np.hstack((xc[1]-xc[-1], xc[2:] - xc[:-2], xc[0]-xc[-2]))
dyc  = np.hstack((yc[1]-yc[-1], yc[2:] - yc[:-2], yc[0]-yc[-2]))
dl = np.sqrt(dxc**2 + dyc**2)
nx, ny = divide(-dyc, dl), divide(dxc, dl)

def trapint(x, y):
    ym = (y[1:] + y[:-1])/2
    xm = (x[1:] + x[:-1])/2
    return xm, np.cumsum(ym * np.diff(x))

def intdS(f, mask=lake):
    f = f[1:-1, 1:-1]
    f = (f*dx).T
    f = (f*dy).T
    return f[lake[1:-1, 1:-1]].sum()

def intdl(f):
    return (f*dl).sum()

for i, (q, t) in enumerate(zip(solutions, times)):
    print(i, t)
    fig, (ax1, ax2) = plt.subplots(ncols=2)
    # Lake
    h, hu, hv, s = q
    eta = s - TOPM["lake_alt"]
    u = divide(hu, h)
    v = divide(hv, h)
    u2 = u**2+v**2
    eta = np.where(h < 1e-10, 0., eta)
    # lake = (s-h < (TOPM["lake_alt"]+2)) & (h > 1e-10)
    ax1.imshow(np.ma.MaskedArray(s, mask=~lake), origin="lower", extent=extent)
    Elake[i] = 1/2*rhow*intdS(g*eta**2 + h*u2)
    Vlak[i] = intdS(h, mask=lake)
    # Avalanche TSUL
    h, hu, hv, s = read_TSUL_contour(i, xc, yc)
    eta = s - TOPM["lake_alt"]
    u = divide(hu, h)
    v = divide(hv, h)
    un = u*nx+v*ny
    u2 = u**2+v**2
    ax2.plot(un, ls="-", marker='o', ms=2, mfc="w", label="TSUL")
    Etsul[i] = 1/2*rhow*intdl(h*( g*eta*un + u2*un ))
    Vtav[i] = intdl(h*un)
    print(f"{intdl(un) = }, {np.isnan(un).any()}")
    # Avalanche AVAC
    h, hu, hv, s = interp_AVAC_contour(t, xc, yc)
    eta = s - TOPM["lake_alt"]
    u = divide(hu, h)
    v = divide(hv, h)
    un = u*nx+v*ny
    u2 = u**2+v**2
    ax2.plot(un, ls="-", marker='o', ms=2, mfc="w", label="AVAC")
    Eavac[i] = 1/2*rhos*intdl(h*( g*eta*un + u2*un ))
    Vavc[i] = intdl(h*un)
    print(f"{intdl(un) = }", end="\n\n")
    ax1.plot(extent[:2], extent[2:], c='k', lw=2)
    ax2.legend()
    plt.colorbar(ax1.scatter(xc, yc, c=u2))
    plt.show(block=False)

tm, Eavac = trapint(times, Eavac)
tm, Etsul = trapint(times, Etsul)

fig, (ax, ax2) = plt.subplots(ncols=2)
ax.plot(tm, Etsul, '-o', label=r"$E_\text{TSUL}$")
ax.plot(tm, Eavac, '-o', label=r"$E_\text{AVAC}$")
ax.plot(times, Elake, '-o', mfc='w', label=r"$\Delta E_\text{lake}$")
ax.legend()
ax2.plot(tm, rhow*Vlak[0]+rhow*trapint(times, Vtav)[1], "-o")
ax2.plot(tm, rhow*Vlak[0]+rhos*trapint(times, Vavc)[1], "-o")
ax2.plot(times, rhow*Vlak, "-.")
plt.show()
