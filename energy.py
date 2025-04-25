from pathlib import Path
from argparse import ArgumentParser
from yaml import safe_load
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.path import Path as mPath
from clawpack.geoclaw import fgout_tools
from clawpack.visclaw.gridtools import grid_output_2d
from clawpack.pyclaw.solution import Solution
try:
    import scienceplots
    plt.style.use('science')
except ImportError:
    plt.style.use("seaborn-v0_8-paper")
# matplotlib.use("pgf")
matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
    'text.usetex': True,
#     'pgf.rcfonts': False,
})
np.seterr(all="raise", under="ignore")


def divide(a, b, fill=0., rtol=1e-05, atol=1e-08):
    c = np.empty(a.shape, dtype=np.float64)
    m = np.isclose(b, 0.)
    c[m] = fill
    m = ~m
    c[m] = a[m] / b[m]
    return c

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("-s", "--save", action="store_true")
    parser.add_argument("-o", "--outdir", type=str, nargs="?", default="_output")
    return parser.parse_args()

def read_clawdata(path, sep="=: ", comments="#", skiprows=0):
    clawdata_trans = dict(T=True, F=False)
    clawdata = dict()
    with open(path) as file:
        lines = [line for line in file.readlines()[skiprows:] if sep in line]
        for line in lines:
            value, key = [e for e in line.split(sep) if e]
            key = key.strip()
            if comments in key:
                key = key[:key.find(comments)].strip()
            value = [v for v in value.split() if v]
            for e, element in enumerate(value):
                try:
                    value[e] = eval(element)
                except Exception:
                    value[e] = clawdata_trans.get(element, element)
            clawdata[key] = value[0] if len(value)==1 else value
    return clawdata

def normal_vectors(x, y):
    dx = np.hstack((x[1]-x[-1], x[2:] - x[:-2], x[0]-x[-2]))/2
    dy = np.hstack((y[1]-y[-1], y[2:] - y[:-2], y[0]-y[-2]))/2
    dl = np.sqrt(dx**2 + dy**2)
    nx = divide(-dy, dl)
    ny = divide(+dx, dl)
    return nx, ny, dl

args = parse_args()
TSUL_outdir = Path(args.outdir).expanduser()

TSULdata = read_clawdata(TSUL_outdir / "fgout_grids.data", sep="#", skiprows=7)
TSULsetprob = read_clawdata(TSUL_outdir / "setprob.data")
AVAC_outdir = Path(TSULsetprob["AVAC_outdir"]).expanduser()
AVACdata = read_clawdata(AVAC_outdir / "fgout_grids.data", sep="#", skiprows=7)

output_formats = [None, "ascii", "binary32", "binary64", "HDF5"]  # TODO implement HDF5

TSULgeodata = read_clawdata(TSUL_outdir / "geoclaw.data")
voellmy = read_clawdata(AVAC_outdir / "voellmy.data")
dry_tolerance = TSULgeodata["dry_tolerance"]
g = TSULgeodata["gravity"]
rhow = TSULgeodata["rho"]
rhos = voellmy["snow_density"]

xmin, ymin = TSULdata["x1, y1"]
xmax, ymax = TSULdata["x2, y2"]
numx, numy = TSULdata["nx,ny"]
x = np.linspace(xmin, xmax, num=numx, endpoint=True)
y = np.linspace(ymin, ymax, num=numy, endpoint=True)
dx = np.diff(x).mean()
dy = np.diff(y).mean()
X, Y = np.meshgrid(x, y)

contour = np.loadtxt("contour.xy").T
print("Creating lake mask...", end=" ", flush=True)
lake = mPath(contour.T).contains_points(np.column_stack((X.flatten(), Y.flatten()))).reshape(X.shape)
print("Done.", flush=True)

xc, yc = contour
inside = (xmin <= xc) & (xc <= xmax) & (ymin <= yc) & (yc <= ymax)
xc = xc[inside]
yc = yc[inside]
nx, ny, dl = normal_vectors(xc, yc)


def trapint(x, y):
    ym = (y[1:] + y[:-1])/2
    xm = (x[1:] + x[:-1])/2
    return xm, np.cumsum(ym * np.diff(x))

def intdA(f, x=None, y=None):
    return ((f*dx).T*dy).sum()

def intdl(f, dl):
    return (f*dl).sum()

def read_lake(i, outdir=TSUL_outdir, file_format=output_formats[TSULdata["output_format"]]):
    fgout = np.fromfile(outdir / f"fgout0001.b{i+1:0>4}", np.float64)
    fgout = fgout.reshape(4, numx, numy, order="F")
    return np.swapaxes(fgout, 1, 2)

def read_contour(i, x, y, outdir=TSUL_outdir, file_format=output_formats[TSULdata["output_format"]]):
    frame_sol = Solution(int(i), path=outdir, file_format=file_format)
    return grid_output_2d(frame_sol, lambda q: q, x, y, levels = "all", return_ma=True)

def avalanche_energy_volume(q, rho, h0=0, nx=nx, ny=ny, dl=dl):
    h, hu, hv, s = q
    hun = hu*nx + hv*ny
    u2 = divide(hu**2 + hv**2, h**2)
    pot = 1/2 * rho * g * (s-TSULsetprob["lake_alt"]) * hun
    kin = 1/2 * rho * u2 * hun
    return intdl(pot + kin, dl=dl), intdl(hun, dl=dl), (h-h0).max(), np.sqrt(divide((hu**2+hv**2), h**2)).max()

def lake_energy_volume_alt(q, rho=rhow, h0=0):
    h, hu, hv, _ = q
    hu2 = divide(hu**2 + hv**2, h)
    pot = 1/2 * rho * g * np.where(lake, (h-h0)**2, 0.)
    kin = 1/2 * rho * np.where(lake, hu2, 0.)
    return intdA(pot + kin), intdA(np.where(lake, h, 0.)), (h-h0)[lake0].max(), np.sqrt(divide((hu**2+hv**2), h**2)).max()

def uniform_grid_interp(x, y, Z, xZ=False, yZ=False):
    """Interpolate values on a line (x, y) on a grid Z(var, xZ, yZ)"""
    if xZ is not None:
        x = (x-xZ[0])/(xZ[-1]-xZ[0])*(Z.shape[2]-1)
    if yZ is not None:
        y = (y-yZ[0])/(yZ[-1]-yZ[0])*(Z.shape[1]-1)
    x = np.clip(0., Z.shape[2], x)
    y = np.clip(0., Z.shape[1], y)
    w = np.clip(x.astype(int), 0, Z.shape[2]-2)
    s = np.clip(y.astype(int), 0, Z.shape[1]-2)
    e = w + 1
    n = s + 1
    return (
        + (x-w)*(y-s)*Z[:, n, e]
        + (x-w)*(n-y)*Z[:, s, e]
        + (e-x)*(y-s)*Z[:, n, w]
        + (e-x)*(n-y)*Z[:, s, w]
    )


Elake = np.zeros(TSULdata["nout"], np.float64)
Etsul = np.zeros(TSULdata["nout"], np.float64)
Vtsul = np.zeros(TSULdata["nout"], np.float64)
Vlake = np.zeros(TSULdata["nout"], np.float64)
Eavac = np.zeros(AVACdata["nout"], np.float64)
Vavac = np.zeros(AVACdata["nout"], np.float64)
hlake = np.zeros(TSULdata["nout"], np.float64)
havac = np.zeros(AVACdata["nout"], np.float64)
htsul = np.zeros(TSULdata["nout"], np.float64)
vlake = np.zeros(TSULdata["nout"], np.float64)
vavac = np.zeros(AVACdata["nout"], np.float64)
vtsul = np.zeros(TSULdata["nout"], np.float64)

AVACtimes = []
for i in range(AVACdata["nout"]):
    t = read_clawdata(AVAC_outdir/f"fgout0001.t{i+1:0>4}", sep=" ")["time"]
    AVACtimes.append(t)
    print(f"AVAC: {t = :.1f} ", end=f"({(i+1)/AVACdata['nout']:.1%})... \r", flush=True)
    Eavac[i], Vavac[i], havac[i], vavac[i] = avalanche_energy_volume(read_contour(i, xc, yc, outdir=AVAC_outdir), rhos)
print(flush=True)
t_start = AVACtimes[max(0, (havac > 0).argmax()-1)]

q0 = read_lake(i)
lake0 = (q0[0] > dry_tolerance) & lake
# qi0 = uniform_grid_interp(xc, yc, q0, x, y)

TSULtimes = []
for i in range(TSULdata["nout"]):
    t = read_clawdata(TSUL_outdir/f"fgout0001.t{i+1:0>4}", sep=" ")["time"]
    TSULtimes.append(t)
    print(f"TSUL: {t = :.1f} ", end=f"({(i+1)/TSULdata['nout']:.1%})... \r", flush=True)
    q = read_lake(i)
    Elake[i], Vlake[i], hlake[i], vlake[i] = lake_energy_volume_alt(q, rhow, q0[0])
    Etsul[i], Vtsul[i], htsul[i], vtsul[i] = avalanche_energy_volume(uniform_grid_interp(xc, yc, q, x, y), rhow)
print(flush=True)
t_start = AVACtimes[max(0, (havac > 0).argmax()-1)]

AVACtimes = np.array(AVACtimes)
TSULtimes = np.array(TSULtimes)

ta, Eavac = trapint(AVACtimes, Eavac)
tm, Etsul = trapint(TSULtimes, Etsul)
Eavac += Elake[0]
Etsul += Elake[0]

fig, ((ax, ax3), (ax2, ax4)) = plt.subplots(ncols=2, nrows=2, layout="tight", figsize=(8, 6))
fig.suptitle(TSUL_outdir)

lavac, = ax.plot(ta, Eavac, ':', label="AVAC")
ltsul, = ax.plot(tm, Etsul,  '-.', label="TSUL")
llake, = ax.plot(TSULtimes, Elake, '-', label=r"$\Delta \mathcal{E_L}$")
ax.set_xlabel("$t$ [s]")
ax.set_ylabel(r"$\mathcal{E}$ [J]")
ax.legend()

Elake = (Elake[1:] + Elake[:-1])/2
Vlake = (Vlake[1:] + Vlake[:-1])/2

ax2.plot(ta, rhos*trapint(AVACtimes, Vavac)[1], ls=lavac.get_linestyle(), c=lavac.get_color(), label="AVAC")
ax2.plot(tm, rhow*trapint(TSULtimes, Vtsul)[1], ls=ltsul.get_linestyle(), c=ltsul.get_color(), label="TSUL")
ax2.plot(tm, rhow*(Vlake-Vlake[0]), ls=llake.get_linestyle(), c=llake.get_color(), label=r"$\Delta M_\mathcal{L}$")
ax2.legend()
ax2.set_xlabel("$t$ [s]")
ax2.set_ylabel(r"$M$ [kg]")
ax2.sharex(ax)

ax4.plot(AVACtimes, havac, ls=lavac.get_linestyle(), c=lavac.get_color(), label="AVAC")
ax4.plot(TSULtimes, htsul, ls=ltsul.get_linestyle(), c=ltsul.get_color(), label="TSUL")
ax4.plot(TSULtimes, hlake, ls=llake.get_linestyle(), c=llake.get_color(), label=r"$\mathcal{L}$")
# ax4.plot(TSULtimes, hdam, alpha=1.0, color="k", label=r"$\mathcal{B}$")
ax4.set_xlabel("$t$ [s]")
ax4.set_ylabel(r"$h_{max}$ [m]")
ax4.legend()

# ax5 = ax4.twinx()
ax5 = ax3
ax5.plot(AVACtimes, vavac, ls=lavac.get_linestyle(), color=lavac.get_color(), label="AVAC")
ax5.plot(TSULtimes, vtsul, ls=ltsul.get_linestyle(), color=ltsul.get_color(), label="TSUL")
ax5.plot(TSULtimes, vlake, ls=llake.get_linestyle(), color=llake.get_color(), label=r"$\mathcal{L}$")
# ax5.plot(TSULtimes, vdam, alpha=1.0, color="k", label=r"$\mathcal{B}$")
ax5.set_ylabel(r"$|u|$ [m/s]")
ax5.sharex(ax)
ax5.legend()

plt.show()

