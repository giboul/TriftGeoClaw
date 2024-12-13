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

def divide(a, b, fill=0., rtol=1e-05, atol=1e-08):
    c = np.empty(a.shape)
    m = np.isclose(b, 0.)
    c[m] = fill
    m = ~m
    c[m] = a[m] / b[m]
    return c

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
TSULgeodata = read_clawdata(projdir / "TSUL" / outdir / "geoclaw.data")
voellmy = read_clawdata(projdir / "AVAC" / outdir / "voellmy.data")
dry_tolerance = TSULgeodata["dry_tolerance"]
g = TSULgeodata["gravity"]
rhow = TSULgeodata["rho"]
rhos = voellmy["snow_density"]
xmin, ymin = TSULdata["lower"]
xmax, ymax = TSULdata["upper"]
numx, numy = TSULdata["num_cells"]
times = np.loadtxt(projdir / "TSUL" / outdir / "times.txt")

x = np.linspace(xmin, xmax, numx)
y = np.linspace(ymin, ymax, numy)
X, Y = np.meshgrid(x, y)
extent = x.min(), x.max(), y.min(), y.max()
dx = (xmax-xmin)/(numx-1)
dy = (ymax-ymin)/(numy-1)

# contour1 = np.loadtxt(projdir / "TOPM" / "contour1.xy").T
contour2 = np.loadtxt(projdir / "TOPM" / "contour2.xy").T
lake = mPath(contour2.T).contains_points(np.column_stack((X.flatten(), Y.flatten()))).reshape(X.shape)
xc, yc = contour2
inside = (extent[0] <= xc) & (xc <= extent[1]) & (extent[2] <= yc) & (yc <= extent[3])
xc = xc[inside]
yc = yc[inside]
dxc  = np.hstack((xc[1]-xc[-1], xc[2:] - xc[:-2], xc[0]-xc[-2]))/2
dyc  = np.hstack((yc[1]-yc[-1], yc[2:] - yc[:-2], yc[0]-yc[-2]))/2
dl = np.sqrt(dxc**2 + dyc**2)
nx = divide(-dyc, dl)
ny = divide(+dxc, dl)

Elake = np.zeros(times.size, np.float64)
Etsul = np.zeros(times.size, np.float64)
Eavac = np.zeros(times.size, np.float64)
Vavac = np.zeros(times.size, np.float64)
Vtsul = np.zeros(times.size, np.float64)
Vlake = np.zeros(times.size, np.float64)

def trapint(x, y):
    ym = (y[1:] + y[:-1])/2
    xm = (x[1:] + x[:-1])/2
    return xm, np.cumsum(ym * np.diff(x))

def intdS(f, x=None, y=None):
    return ((f*dx).T*dy).sum()

def intdl(f, x=None, y=None):
    return (f*dl).sum()

def read_lake(i, outdir=projdir/"TSUL"/outdir, file_format=TSUL["out_format"]):
    frame = Solution(int(i), path=outdir, file_format=file_format)
    return grid_output_2d(frame, lambda q: q, X, Y, levels='all', return_ma=True)

def read_contour(i, x, y, outdir=projdir/"TSUL"/outdir, file_format=TSUL["out_format"]):
    frame_sol = Solution(int(i), path=outdir, file_format=file_format)
    return grid_output_2d(frame_sol, lambda q: q, x, y, levels = "all", return_ma=True)

def interp_contour(t, times, x, y, outdir=projdir/"AVAC"/outdir, file_format=AVAC["out_format"]):
    it = np.clip((times <= t).sum()-1, 0, times.size-2)
    q1 = read_contour(int(it), x, y, outdir, file_format=file_format)
    q2 = read_contour(int(it+1), x, y, outdir, file_format=file_format)
    qi = q1 + (q2 - q1) * (t - times[it])/(times[it+1] - times[it])
    return qi

def lake_energy_volume(q, rho):
    h, hu, hv, s = q
    s = s - TOPM["lake_alt"]
    hu2 = divide(hu**2 + hv**2, h)
    pot = 1/2 * rhow * g * np.where(lake, (s-h)**2, 0.)
    kin = 1/2 * rhow * hu2
    print(f"{pot.min(), pot.max() = }")
    return intdS(pot + kin), intdS(np.where(lake, h, 0.))

def avalanche_energy_volume(q, rho):
    h, hu, hv, s = q
    s = s - TOPM["lake_alt"]
    hun = hu*nx + hv*ny
    un = divide(hun, h)
    u2 = divide(hu**2 + hv**2, h**2)
    pot = 1/2 * rho * g * (s**2-(s-h)**2) * un
    kin = 1/2 * rho * u2 * hun
    return intdl(pot + kin), intdl(hun)

for i, t in enumerate(times):
    Eavac[i], Vavac[i] = avalanche_energy_volume(interp_contour(t, AVACtimes, xc, yc), rhos)
    Etsul[i], Vtsul[i] = avalanche_energy_volume(read_contour(i, xc, yc), rhow)
    Elake[i], Vlake[i] = lake_energy_volume(read_lake(i), rhow)

tm, Eavac = trapint(times, Eavac)
tm, Etsul = trapint(times, Etsul)
print(Elake[0])
Eavac += Elake[0]
Etsul += Elake[0]

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
