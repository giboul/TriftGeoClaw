from pathlib import Path
from argparse import ArgumentParser
from yaml import safe_load
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.path import Path as mPath
from clawpack.visclaw.gridtools import grid_output_2d
from clawpack.pyclaw.solution import Solution
from skimage.morphology import isotropic_erosion
plt.style.use("bmh")


projdir = Path(__file__).parent
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    AVAC = config["AVAC"]
    TOPM = config["TOPM"]
    TSUL = config["TSUL"]

def divide(a, b, fill=0., rtol=1e-05, atol=1e-08):
    c = np.empty(a.shape, dtype=np.float64)
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
TSULdir = projdir / "TSUL" / outdir
AVACdir = projdir / "AVAC" / outdir

AVACdata = read_clawdata(AVACdir / "claw.data")
AVACtimes = np.linspace(AVACdata["t0"], AVACdata["tfinal"], AVACdata["num_output_times"], endpoint=True)

TSULdata = read_clawdata(TSULdir / "claw.data")
TSULtimes = np.linspace(TSULdata["t0"], TSULdata["tfinal"], TSULdata["num_output_times"], endpoint=True)
TSULgeodata = read_clawdata(TSULdir / "geoclaw.data")
voellmy = read_clawdata(AVACdir / "voellmy.data")
dry_tolerance = TSULgeodata["dry_tolerance"]
g = TSULgeodata["gravity"]
rhow = TSULgeodata["rho"]
rhos = voellmy["snow_density"]

xmin, ymin = TSULdata["lower"]
xmax, ymax = TSULdata["upper"]
numx, numy = TSULdata["num_cells"]

# times = np.loadtxt(TSULdir / "times.txt")
# files = list(sorted(TSULdir.glob("fort.t*")))
# times = np.empty(len(files), dtype=np.float64)
# for i, path in enumerate(files):
#     with open(path) as file:
#         line = file.readline()
#         times[i] = float(line[:len(line)-line[::-1].find(" ")])

x = np.linspace(xmin, xmax, numx)
y = np.linspace(ymin, ymax, numy)
X, Y = np.meshgrid(x, y)
extent = x.min(), x.max(), y.min(), y.max()
dx = (xmax-xmin)/(numx-1)
dy = (ymax-ymin)/(numy-1)

# contour1 = np.loadtxt(projdir / "TOPM" / "contour1.xy").T
contour2 = np.loadtxt(projdir / "TOPM" / "contour2.xy").T
lake = mPath(contour2.T).contains_points(np.column_stack((X.flatten(), Y.flatten()))).reshape(X.shape)
# lake = isotropic_erosion(lake, 1)
xc, yc = contour2
inside = (extent[0] <= xc) & (xc <= extent[1]) & (extent[2] <= yc) & (yc <= extent[3])
xc = xc[inside]
yc = yc[inside]
dxc  = np.hstack((xc[1]-xc[-1], xc[2:] - xc[:-2], xc[0]-xc[-2]))/2
dyc  = np.hstack((yc[1]-yc[-1], yc[2:] - yc[:-2], yc[0]-yc[-2]))/2
dl = np.sqrt(dxc**2 + dyc**2)
nx = divide(-dyc, dl)
ny = divide(+dxc, dl)

Elake = np.zeros(TSULtimes.size, np.float64)
Etsul = np.zeros(TSULtimes.size, np.float64)
Vtsul = np.zeros(TSULtimes.size, np.float64)
Vlake = np.zeros(TSULtimes.size, np.float64)
Eavac = np.zeros(AVACtimes.size, np.float64)
Vavac = np.zeros(AVACtimes.size, np.float64)

def trapint(x, y):
    ym = (y[1:] + y[:-1])/2
    xm = (x[1:] + x[:-1])/2
    return xm, np.cumsum(ym * np.diff(x))

def intdA(f, x=None, y=None):
    return ((f*dx).T*dy).sum()

def intdl(f, x=None, y=None):
    return (f*dl).sum()

def read_lake(i, outdir=TSULdir, file_format=TSUL["out_format"]):
    frame = Solution(int(i), path=outdir, file_format=file_format)
    return grid_output_2d(frame, lambda q: q, X, Y, levels='all', return_ma=True)

def read_contour(i, x, y, outdir=TSULdir, file_format=TSUL["out_format"]):
    frame_sol = Solution(int(i), path=outdir, file_format=file_format)
    return grid_output_2d(frame_sol, lambda q: q, x, y, levels = "all", return_ma=True)

# def interp_contour(t, times, x, y, outdir=AVACdir, file_format=AVAC["out_format"]):
#     it = np.clip((times <= t).sum()-1, 0, times.size-2)
#     q1 = read_contour(int(it), x, y, outdir, file_format=file_format)
#     q2 = read_contour(int(it+1), x, y, outdir, file_format=file_format)
#     qi = q1 + (q2 - q1) * (t - times[it])/(times[it+1] - times[it])
#     return qi

h0, u0, v0, s0 = read_lake(0)

def lake_energy_volume(q, rho):
    h, hu, hv, s = q
    hu2 = divide(hu**2 + hv**2, h)
    pot = 1/2 * rhow * g * np.where(lake, (h-h0)**2, 0.)
    kin = 1/2 * rhow * np.where(lake, hu2, 0.)
    # f, (a1, a2, a3) = plt.subplots(ncols=3)
    # a1.imshow(np.where(lake, h, np.nan), origin="lower", cmap=plt.cm.RdBu)
    # a2.imshow(np.where(lake, pot, np.nan), origin="lower", cmap=plt.cm.RdBu)
    # i = a3.imshow(np.where(lake, s-s0, np.nan), origin="lower", cmap=plt.cm.RdBu)
    # plt.colorbar(i, ax=a3)
    # i.set_clim(-np.abs(s-s0)[lake].max(), np.abs(s-s0)[wet_lake].max())
    # plt.show()
    return intdA(pot + kin), intdA(np.where(lake, h, 0.))

def avalanche_energy_volume(q, rho):
    h, hu, hv, s = q
    hun = hu*nx + hv*ny
    u2 = divide(hu**2 + hv**2, h**2)
    pot = 1/2 * rho * g * (s-TOPM["lake_alt"]) * hun
    kin = 1/2 * rho * u2 * hun
    return intdl(pot + kin), intdl(hun)

for i, t in enumerate(AVACtimes):
    Eavac[i], Vavac[i] = avalanche_energy_volume(read_contour(i, xc, yc, outdir=AVACdir), rhos)
for i, t in enumerate(TSULtimes):
    Etsul[i], Vtsul[i] = avalanche_energy_volume(read_contour(i, xc, yc), rhow)
    Elake[i], Vlake[i] = lake_energy_volume(read_lake(i), rhow)

ta, Eavac = trapint(AVACtimes, Eavac)
tm, Etsul = trapint(TSULtimes, Etsul)
Eavac += Elake[0]
Etsul += Elake[0]

fig, (ax, ax2, ax3) = plt.subplots(ncols=3, layout="tight")
lavac, = ax.plot(ta, Eavac, '-.o', label=r"$\mathcal{E_A}$ (AVAC)", ms=3, mfc='w')
ltsul, = ax.plot(tm, Etsul,  '-o', label=r"$\mathcal{E_A}$ (TSUL)", ms=3, mfc='w')
llake, = ax.plot(TSULtimes, Elake, '-o', label=r"$\mathcal{\Delta E_L}$ (TSUL)", ms=3, mfc='w')

Elake = (Elake[1:] + Elake[:-1])/2
Vlake = (Vlake[1:] + Vlake[:-1])/2

ax.legend()
ax.set_xlabel("$t$ [s]")
ax.set_ylabel(r"$\mathcal{E}$ [J]")
ax2.plot(ta, rhos*trapint(AVACtimes, Vavac)[1], ls=lavac.get_linestyle(), c=lavac.get_color(), label=r"$V_\mathcal{A}$ (AVAC)")
ax2.plot(tm, rhow*trapint(TSULtimes, Vtsul)[1], ls=ltsul.get_linestyle(), c=ltsul.get_color(), label=r"$V_\mathcal{A}$ (TSUL)")
ax2.plot(tm, rhow*(Vlake-Vlake[0]), ls=llake.get_linestyle(), c=llake.get_color(), label=r"$\Delta V_\mathcal{L}$ (TSUL)")
ax2.legend()
ax2.set_xlabel("$t$ [s]")
ax2.set_ylabel(r"$V_\mathcal{L}$ [m$^\mathregular{3}$]")
ax3.plot(np.interp(tm, ta, Eavac), Elake, ls=lavac.get_linestyle(), c=lavac.get_color(), label="AVAC")
ax3.plot(Etsul, Elake, ls=ltsul.get_linestyle(), c=ltsul.get_color(), label="TSUL")
ax3.axline((0, 0), slope=1, ls='-.', c="gray", label=r"$\mathcal{E_A=\Delta E_L}$")
ax3.set_xlabel(r"$\mathcal{E_L}$ [J]")
ax3.set_ylabel(r"$\mathcal{E_A}$ [J]")
ax3.legend()
plt.show()
