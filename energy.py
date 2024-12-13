from pathlib import Path
from argparse import ArgumentParser
from yaml import safe_load
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.path import Path as mPath
from clawpack.visclaw.gridtools import grid_output_2d
from clawpack.pyclaw.solution import Solution


projdir = Path(__file__).parent
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    AVAC = config["AVAC"]
    TOPM = config["TOPM"]
    TSUL = config["TSUL"]

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("avid", type=str, nargs="?", default="")
    return parser.parse_args()

def read_clawdata(path):
    clawdata_trans = dict(T=True, F=False)
    clawdata = dict()
    with open(path) as file:
        lines = [l for l in file.readlines() if "=" in l]
        for line in lines:
            value, key = line.split("=: ")
            key = key[:key.find(" ")]
            if "#" in value:
                value = value[value.find("#")]
            value = [v for v in value.split() if v]
            for e, element in enumerate(value):
                try:
                    value[e] = eval(element)
                except Exception:
                    value[e] = clawdata_trans[element]
            clawdata[key] = value[0] if len(value)==1 else value
    return clawdata

outdir = Path(f"_output{parse_args().avid}")

AVACdata = read_clawdata(projdir / "AVAC" / outdir / "claw.data")
AVACtimes = np.linspace(AVACdata["t0"], AVACdata["tfinal"], AVACdata["num_output_times"], endpoint=True)

TSULdata = read_clawdata(projdir / "TSUL" / outdir / "claw.data")
dry_tolerance = read_clawdata(projdir / "TSUL" / outdir / "geoclaw.data")["dry_tolerance"]
xmin, ymin = TSULdata["lower"]
xmax, ymax = TSULdata["upper"]
numx, numy = TSULdata["num_cells"]
times = np.loadtxt(projdir/"TSUL"/outdir/"times.txt")

x = np.linspace(xmin, xmax, 100)
y = np.linspace(ymin, ymax, 100)
X, Y = np.meshgrid(x, y)
extent = X.min(), X.max(), Y.min(), Y.max()
# dx = (xmax-xmin)/numx
# dy = (ymax-ymin)/numy
dx = (x[2:] - x[:-2])/2
dy = (y[2:] - y[:-2])/2

contour1 = np.loadtxt(projdir / "TOPM" / "contour1.xy").T
contour2 = np.loadtxt(projdir / "TOPM" / "contour2.xy").T


def read_TSUL_lake(i, outdir=projdir/"TSUL"/outdir):
    frame = Solution(int(i), path=outdir, file_format=TSUL["out_format"])
    return grid_output_2d(frame, lambda x: x, X, Y, levels='all', return_ma=True)


def read_TSUL_contour(i, x, y, outdir=projdir/"TSUL"/outdir):
    frame_sol = Solution(int(i), path=outdir, file_format=TSUL["out_format"])
    return grid_output_2d(frame_sol, lambda q: q, x, y, levels = "all", return_ma=True)


def interp_AVAC_contour(t, x, y, outdir=projdir/"AVAC"/outdir):
    it = np.clip((AVACtimes <= t).sum()-1, 0, AVACtimes.size-2)
    frame_sol = Solution(int(it), path=outdir, file_format=AVAC["out_format"])
    q1 = grid_output_2d(frame_sol, lambda q: q, x, y, levels = "all", return_ma=True)
    frame_sol = Solution(int(it+1), path=outdir, file_format=AVAC["out_format"])
    q2 = grid_output_2d(frame_sol, lambda q: q, x, y, levels = "all", return_ma=True)
    qi = q1 + (q2 - q1) * (t - AVACtimes[it])/(AVACtimes[it+1] - AVACtimes[it])
    return qi


def divide(a, b, fill=0., rtol=1e-05, atol=1e-08):
    c = np.empty(a.shape)
    m = np.isclose(b, 0.)
    c[m] = fill
    m = ~m
    c[m] = a[m] / b[m]
    return c


g = 9.81
rhow = 1000
rhos = 300
Elake = np.zeros(times.size, np.float64)
Etsul = np.zeros(times.size, np.float64)
Eavac = np.zeros(times.size, np.float64)
Vavac = np.zeros(times.size, np.float64)
Vtsul = np.zeros(times.size, np.float64)
Vlake = np.zeros(times.size, np.float64)

h0, _, _, eta0 = read_TSUL_lake(0)
z = eta0 - h0
lake = h0 >= dry_tolerance
lake = mPath(contour2.T).contains_points(np.column_stack((X.flatten(), Y.flatten()))).reshape(X.shape)
xc, yc = contour2
inside = (extent[0] <= xc) & (xc <= extent[1]) & (extent[2] <= yc) & (yc <= extent[3])
xc = xc[inside]
yc = yc[inside]
dxc  = np.hstack((xc[1]-xc[-1], xc[2:] - xc[:-2], xc[0]-xc[-2]))/2
dyc  = np.hstack((yc[1]-yc[-1], yc[2:] - yc[:-2], yc[0]-yc[-2]))/2
dl = np.sqrt(dxc**2 + dyc**2)
nx = np.zeros_like(dyc)
ny = np.zeros_like(dyc)
np.divide(-dyc, dl, out=nx, where=~np.isclose(dl, 0., atol=1e-10))
np.divide(+dxc, dl, out=ny, where=~np.isclose(dl, 0., atol=1e-10))

def trapint(x, y):
    ym = (y[1:] + y[:-1])/2
    xm = (x[1:] + x[:-1])/2
    return xm, np.cumsum(ym * np.diff(x))

def intdS(f, x=None, y=None):
    if x is not None:
        d_x = (x[2:] - x[:-2])/2
    else:
        d_x = dx
    if y is not None:
        d_y = (y[2:] - y[:-2])/2
    else:
        d_y = dy
    f = f[1:-1, 1:-1]
    f = ((f*d_x).T*d_y).T
    return f.sum()

def intdl(f, x=None, y=None):
    if x is not None and y is not None:
        dx  = np.hstack((x[1]-x[-1], x[2:] - x[:-2], x[0]-x[-2]))/2
        dy  = np.hstack((y[1]-y[-1], y[2:] - y[:-2], y[0]-y[-2]))/2
        d_l = np.sqrt(dx**2 + dy**2)
    else:
        d_l = dl
    return (f*d_l).sum()

def lake_energy(i):
    h, hu, hv, s = read_TSUL_lake(i)
    wet = h > dry_tolerance
    ax1.imshow(np.where(wet&lake, s, np.nan), origin="lower", extent=extent)
    s[~wet|~lake] = TOPM["lake_alt"]
    hu2 = np.zeros_like(hu, dtype=np.float64)
    np.divide(hu**2 + hv**2, h, out=hu2, where=wet)
    Elake[i] = 1/2 * rhow * intdS(g*(s-TOPM["lake_alt"])**2 + hu2)
    Vlake[i] = intdS(np.where(lake, h, 0.))
    print(f"{Vlake[i] = :.3e}")
    ax3.imshow(np.where(wet, hu**2+hv**2, np.nan), origin="lower", extent=extent)
    h, hu, hv, s = read_TSUL_lake(i, outdir=projdir/"AVAC"/outdir)

def water_avalanche_energy(i):
    h, hu, hv, s = read_TSUL_contour(i, xc, yc)
    eta = s - TOPM["lake_alt"]
    eta = np.where(h < dry_tolerance, 0., eta)
    hun = hu*nx + hv*ny
    _hunorm = np.sqrt(hu**2+hv**2).max() or 1.
    for _x, _y, _nx, _ny, _hu, _hv, _hun in zip(xc, yc, nx, ny, hu, hv, hun):
        ax1.arrow(_x, _y, _hu/_hunorm*1e3, _hv/_hunorm*1e3, head_width=10)
        ax3.arrow(_x, _y, _nx*1e2,_ny*1e2, head_width=10)
    ax3.plot(xc, yc)
    u2 = divide(hu**2 + hv**2, h**2)
    ax2.plot(np.cumsum(dl), hun/np.abs(hun).max(), '-o', mfc='w', label="TSUL: $hu_n$")
    ax2.set_title(intdl(hun))
    Etsul[i] = rhow * intdl((g*eta + 1/2*u2)*hun)
    Vtsul[i] = intdl(hun)
    ax1.scatter(xc, yc, c=hun, s=10, vmin=-np.std(hun), vmax=np.std(hun), cmap=plt.cm.RdBu)

def snow_avalanche_energy(t):
    h, hu, hv, s = interp_AVAC_contour(t, xc, yc)
    eta = s - TOPM["lake_alt"]
    u = divide(hu, h)
    v = divide(hv, h)
    un = u*nx + v*ny
    u2 = u**2 + v**2
    hun = hu*nx + hv*ny
    u2 = divide(hu**2 + hv**2, h**2)
    Eavac[i] = 1/2 * rhos * intdl(hun*( g*eta + 1/2*u2 ))
    Vavac[i] = intdl(h*un)
    ax3.scatter(xc, yc, c=u2, s=1)
    ax2.plot(np.cumsum(dl), hun/np.abs(hun).max(), '-o', mfc='w', label="AVAC: $hu_n$")

for i, t in enumerate(times):
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3)
    ax3.sharex(ax1)
    ax3.sharey(ax1)
    lake_energy(i)
    snow_avalanche_energy(t)
    water_avalanche_energy(i)
    ax2.legend()
    plt.show(block=False)

tm, Eavac = trapint(times, Eavac)
tm, Etsul = trapint(times, Etsul)
Elake = Elake - Elake[0]

fig, (ax, ax2, ax3) = plt.subplots(ncols=3)
ax.plot(tm, Etsul, '-o', label=r"$E_\text{TSUL}$")
ax.plot(tm, Eavac, '-o', label=r"$E_\text{AVAC}$")
ax.plot(times, Elake, '-o', mfc='w', label=r"$\Delta E_\text{lake}$")
ax.legend()
ax2.plot(tm, rhow*Vlake[0]+rhow*trapint(times, Vtsul)[1], "-o")
ax2.plot(tm, rhow*Vlake[0]+rhos*trapint(times, Vavac)[1], "-o")
ax2.plot(times, rhow*Vlake, "-.")
ax2.yaxis.tick_right()
ax3.axline((0, 0), slope=1, ls='-.', c="gray")
ax3.plot(Etsul, Elake[1:])
ax3.plot(Eavac, Elake[1:])
plt.show()
