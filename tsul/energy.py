from pathlib import Path
from argparse import ArgumentParser
from yaml import safe_load
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.path import Path as mPath
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


projdir = Path(__file__).parent.parent
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
    parser.add_argument("-s", "--save", action="store_true")
    return parser.parse_args()

def read_clawdata(path, sep="=: ", comments="#"):
    clawdata_trans = dict(T=True, F=False)
    clawdata = dict()
    with open(path) as file:
        lines = [line for line in file.readlines() if sep in line]
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

args = parse_args()
outdir = Path(f"_output")
TSULdir = projdir / "tsul" / outdir
AVACdir = projdir / "avac" / outdir

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

read_clawdata(AVACdir/"fgout0001.t0001", sep="  ")

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

x = np.linspace(xmin, xmax, 4*numx, endpoint=True)
y = np.linspace(ymin, ymax, 4*numy, endpoint=True)
extent = x.min(), x.max(), y.min(), y.max()
dx = np.diff(x).mean()
dy = np.diff(y).mean()
X, Y = np.meshgrid(x, y)

contour = np.loadtxt(projdir / "tsul" / "contour.xy").T
lake = mPath(contour.T).contains_points(np.column_stack((X.flatten(), Y.flatten()))).reshape(X.shape)
# xd, yd = np.loadtxt(projdir / "topm" / "dam.xy").T

xc, yc = contour
inside = (extent[0] <= xc) & (xc <= extent[1]) & (extent[2] <= yc) & (yc <= extent[3])
xc = xc[inside]
yc = yc[inside]

def normal_vectors(x, y):
    dx = np.hstack((x[1]-x[-1], x[2:] - x[:-2], x[0]-x[-2]))/2
    dy = np.hstack((y[1]-y[-1], y[2:] - y[:-2], y[0]-y[-2]))/2
    dl = np.sqrt(dx**2 + dy**2)
    nx = divide(-dy, dl)
    ny = divide(+dx, dl)
    return nx, ny, dl

nx, ny, dl = normal_vectors(xc, yc)
# nxd, nyd, dld = normal_vectors(xd, yd)    

# xd, yd = np.loadtxt(projdir / "TOPM" / "dam.xy").T
# xd = xd + 20
# yd = yd - 20

Elake = np.zeros(TSULtimes.size, np.float64)
Etsul = np.zeros(TSULtimes.size, np.float64)
Vtsul = np.zeros(TSULtimes.size, np.float64)
Vlake = np.zeros(TSULtimes.size, np.float64)
Eavac = np.zeros(AVACtimes.size, np.float64)
Vavac = np.zeros(AVACtimes.size, np.float64)
hlake = np.zeros(TSULtimes.size, np.float64)
havac = np.zeros(AVACtimes.size, np.float64)
htsul = np.zeros(TSULtimes.size, np.float64)
vlake = np.zeros(TSULtimes.size, np.float64)
vavac = np.zeros(AVACtimes.size, np.float64)
vtsul = np.zeros(TSULtimes.size, np.float64)
hdam = np.zeros(TSULtimes.size, np.float64)
# Vdam = np.zeros(TSULtimes.size, np.float64)
# vdam = np.zeros(TSULtimes.size, np.float64)
# Edam = np.zeros(TSULtimes.size, np.float64)
# hdama = np.zeros(AVACtimes.size, np.float64)


def trapint(x, y):
    ym = (y[1:] + y[:-1])/2
    xm = (x[1:] + x[:-1])/2
    return xm, np.cumsum(ym * np.diff(x))

def intdA(f, x=None, y=None):
    return ((f*dx).T*dy).sum()

def intdl(f, dl=dl):
    return (f*dl).sum()

def read_lake(i, outdir=TSULdir, file_format=TSUL["out_format"]):
    frame = Solution(int(i), path=outdir, file_format=file_format)
    return grid_output_2d(frame, lambda q: q, X, Y, levels='all', return_ma=True)

def read_contour(i, x, y, outdir=TSULdir, file_format=TSUL["out_format"]):
    frame_sol = Solution(int(i), path=outdir, file_format=file_format)
    return grid_output_2d(frame_sol, lambda q: q, x, y, levels = "all", return_ma=True)

def avalanche_energy_volume(q, rho, h0=0, nx=nx, ny=ny, dl=dl):
    h, hu, hv, s = q
    hun = hu*nx + hv*ny
    u2 = divide(hu**2 + hv**2, h**2)
    pot = 1/2 * rho * g * (s-TOPM["lake_alt"]) * hun
    kin = 1/2 * rho * u2 * hun
    return intdl(pot + kin, dl=dl), intdl(hun, dl=dl), (h-h0).max(), np.sqrt(divide((hu**2+hv**2), h**2)).max()

for i, t in enumerate(AVACtimes):
    print(f"AVAC: {t = :.1f} ", end=f"({(i+1)/AVACtimes.size:.1%})... \r")
    Eavac[i], Vavac[i], havac[i], vavac[i] = avalanche_energy_volume(read_contour(i, xc, yc, outdir=AVACdir), rhos)
print()
t_start = AVACtimes[max(0, (havac > 0).argmax()-1)]

def lake_energy_volume_alt(q, rho=rhow):
    h, hu, hv, _ = q
    hu2 = divide(hu**2 + hv**2, h)
    pot = 1/2 * rho * g * np.where(lake, (h-h0)**2, 0.)
    kin = 1/2 * rho * np.where(lake, hu2, 0.)
    return intdA(pot + kin), intdA(np.where(lake, h, 0.)), (h-h0)[lake0].max(), np.sqrt(divide((hu**2+hv**2), h**2)).max()

if __name__ == "__main__":

    h0 = read_lake(0)[0]
    ht0 = read_contour(0, xc, yc)[0]
    # hd0 = read_contour(0, xd, yd)[0]
    lake0 = (h0 > 1e-5) | lake

    for i, t in enumerate(TSULtimes):
        print(f"TSUL: {t = :.1f} ", end=f"({(i+1)/TSULtimes.size:.1%})... \r")
        Etsul[i], Vtsul[i], htsul[i], vtsul[i] = avalanche_energy_volume(read_contour(i, xc, yc), rhow, h0=ht0)
        Elake[i], Vlake[i], hlake[i], vlake[i] = lake_energy_volume_alt(read_lake(i), rhow)
        # Edam[i],  Vdam[i],  hdam[i],  vdam[i]  = avalanche_energy_volume(read_contour(i, xd, yd, outdir=TSULdir), rhow, h0=hd0, nx=nxd, ny=nyd, dl=dld)
    print()

    ta, Eavac = trapint(AVACtimes, Eavac)
    tm, Etsul = trapint(TSULtimes, Etsul)
    Eavac += Elake[0]
    Etsul += Elake[0]

    fig, ((ax, ax3), (ax2, ax4)) = plt.subplots(ncols=2, nrows=2, layout="tight", figsize=(8, 6))
    fig.suptitle(f"{AVACdir.relative_to(AVACdir.parents[1])} versus {TSULdir.relative_to(TSULdir.parents[1])}")

    lavac, = ax.plot(ta, Eavac, ':', label="AVAC")
    ltsul, = ax.plot(tm, Etsul,  '-.', label="TSUL")
    llake, = ax.plot(TSULtimes, Elake, '-', label=r"$\Delta \mathcal{E_L}$")
    ax.set_xlabel("$t$ [s]")
    ax.set_ylabel(r"$\mathcal{E}$ [J]")
    ax.legend()

    Elake = (Elake[1:] + Elake[:-1])/2
    Vlake = (Vlake[1:] + Vlake[:-1])/2

    # ax3.plot(np.interp(tm, ta, Eavac), Elake, ls=lavac.get_linestyle(), c=lavac.get_color(), label="AVAC")
    # ax3.plot(Etsul, Elake, ls=ltsul.get_linestyle(), c=ltsul.get_color(), label="TSUL")
    # ax3.axline((0, 0), slope=1, ls=llake.get_linestyle(), c=llake.get_color(), label=r"$\mathcal{E_A}=\Delta \mathcal{E_L}$")
    # ax3.set_ylabel(r"$\mathcal{E_L}$ [J]")
    # ax3.set_xlabel(r"$\mathcal{E_A}$ [J]")
    # ax3.legend()

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

    if args.save:
        fig.savefig(projdir / "figures" / f"energy.pdf", bbox_inches="tight")
        # fig.savefig(projdir / "figures" / f"energy.pgf", bbox_inches="tight")
    else:
        plt.show()

    # plt.plot(AVACtimes, hdama)
    # plt.plot(TSULtimes, hdamt)
    # plt.show()

    with open("max.log", "a") as file:
        txt = rf"{vtsul.max(): >5.1f} & {havac.max(): >5.1f} & {Vavac.max():.3e} & {hlake.max(): >5.1f}\\""\n"
        file.write(txt)
