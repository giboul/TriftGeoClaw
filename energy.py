from pathlib import Path
from argparse import ArgumentParser
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.path import Path as mPath
from topo_utils import find_contour
from clawpack.clawutil.data import ClawData
from clawpack.visclaw.gridtools import grid_output_2d
from clawpack.pyclaw.solution import Solution
from clawpack.geoclaw import fgout_tools


def divide(a, b, fill=0., rtol=1e-05, atol=1e-08):
    c = np.empty(a.shape, dtype=np.float64)
    m = np.isclose(b, 0., rtol, atol)
    c[m] = fill
    m = ~m
    c[m] = a[m] / b[m]
    return c

def normal_vectors(x, y):
    dx = np.hstack((x[1]-x[-1], x[2:] - x[:-2], x[0]-x[-2]))/2
    dy = np.hstack((y[1]-y[-1], y[2:] - y[:-2], y[0]-y[-2]))/2
    dl = np.sqrt(dx**2 + dy**2)
    nx = divide(-dy, dl)
    ny = divide(+dx, dl)
    return nx, ny, dl

def trapint(x, y):
    ym = (y[1:] + y[:-1])/2
    xm = (x[1:] + x[:-1])/2
    return xm, np.cumsum(ym * np.diff(x))

def read_contour(i, x, y, outdir, file_format):
    frame_sol = Solution(int(i), path=outdir, file_format=file_format)
    return grid_output_2d(frame_sol, lambda q: q, x, y, levels = "all", return_ma=True)

def avalanche_energy_volume(q, rho, ref_alt, nx, ny, dl, h0=0, g=9.81):
    h, hu, hv, s = q
    hun = hu*nx + hv*ny
    u2 = divide(hu**2 + hv**2, h**2)
    pot = 1/2 * rho * g * (s-ref_alt) * hun
    kin = 1/2 * rho * u2 * hun
    return np.nansum((pot + kin)*dl), np.nansum(hun*dl)

def lake_energy_volume(q, rho, dx, dy, h0=0, g=9.81):
    h, hu, hv, _ = q
    hu2 = divide(hu**2 + hv**2, h)
    pot = 1/2 * rho * g * (h-h0)**2
    kin = 1/2 * rho * np.where(h>0, hu2, 0.)
    return np.nansum((pot + kin)*dx*dy), np.nansum(h*dx*dy)

class Data:
    """Read all datafiles in an output directory at once."""
    def __init__(self, outdir: Path):
        datafiles = list(Path(outdir).glob("*.data"))
        for path in datafiles:
            data = ClawData()
            data.read(path, force=True)
            setattr(self, path.stem, data)


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


def main():

    args = parse_args()
    wave_outdir = Path(args.outdir).expanduser()

    wave_data = Data(wave_outdir)
    avac_outdir = Path(wave_data.setprob.AVAC_outdir.strip("'")).expanduser()

    avac_data = Data(avac_outdir)

    wave_fg = fgout_tools.FGoutGrid(args.fgno, wave_outdir)
    wave_fg.read_fgout_grids_data()

    avac_fg = fgout_tools.FGoutGrid(args.fgno, avac_outdir)
    avac_fg.read_fgout_grids_data()

    wave_fgout_init = wave_fg.read_frame(1)
    lake = wave_fgout_init.h > 0
    lake = wave_fgout_init.B > wave_data.setprob.min_alt_avac
    lake_alt = wave_fgout_init.eta[wave_fgout_init.h>0].mean()
    wave_fg.extent = wave_fg.x1, wave_fg.x2, wave_fg.y1, wave_fg.y2
    xc, yc = find_contour(lake, wave_fg.extent, wave_fg.nx, wave_fg.ny).T
    nx, ny, dl = normal_vectors(xc, yc)

    Ea = np.zeros(avac_fg.nout, np.float64)
    Va = np.zeros(avac_fg.nout, np.float64)
    for i, t in enumerate(avac_fg.times):
        q = avac_fg.read_frame(i+1).q
        qc = uniform_grid_interp(xc, yc, q, avac_fg.x, avac_fg.y)
        Ea[i], Va[i] = avalanche_energy_volume(qc, avac_data.geoclaw.rho, lake_alt, nx, ny, dl)
    print(Ea)
    print(Va)

    Ew = np.zeros(wave_fg.nout, np.float64)
    Vw = np.zeros(wave_fg.nout, np.float64)
    El = np.zeros(wave_fg.nout, np.float64)
    Vl = np.zeros(wave_fg.nout, np.float64)
    dx = (wave_fg.x2 - wave_fg.x1) / wave_fg.nx
    dy = (wave_fg.y2 - wave_fg.y1) / wave_fg.ny
    for i, t in enumerate(wave_fg.times):
        q = wave_fg.read_frame(i+1).q
        qc = uniform_grid_interp(xc, yc, q, wave_fg.x, wave_fg.y)
        Ew[i], Vw[i] = avalanche_energy_volume(qc, wave_data.geoclaw.rho, lake_alt, nx, ny, dl)
        El[i], Vl[i] = lake_energy_volume(qc, wave_data.geoclaw.rho, lake_alt, dx, dy)
    print(Ew)
    print(Vw)
    print(El)
    print(Vl)

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("-s", "--save", action="store_true")
    parser.add_argument("-o", "--outdir", type=str, nargs="?", default="_output")
    parser.add_argument("--fgno", type=int, nargs="?", default=1)
    return parser.parse_args()

if __name__ == "__main__":
    main()
exit()

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


TSULdata = ClawData()
TSULdata.read(TSUL_outdir / "fgout_grids.data", sep="#", skiprows=7)
TSULsetprob = ClawData()
TSULsetprob.read(TSUL_outdir / "setprob.data")
AVACdata = ClawData()
AVAC_outdir = Path(TSULsetprob["AVAC_outdir"]).expanduser()
AVACdata.read(AVAC_outdir / "fgout_grids.data", sep="#", skiprows=7)


TSULgeodata = ClawData()
voellmy = ClawData()
TSULgeodata.read(TSUL_outdir / "geoclaw.data")
voellmy.read(AVAC_outdir / "voellmy.data")
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
    AVACtimes.append(t)
    print(f"AVAC: {t = :.1f} ", end=f"({(i+1)/AVACdata['nout']:.1%})... \r", flush=True)
    Eavac[i], Vavac[i], havac[i], vavac[i] = avalanche_energy_volume(read_contour(i, xc, yc, outdir=AVAC_outdir), rhos)
print(flush=True)
t_start = AVACtimes[max(0, (havac > 0).argmax()-1)]

q0 = read_lake(i)
lake0 = (q0[0] > dry_tolerance) & lake

TSULtimes = []
for i in range(TSULdata["nout"]):
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

