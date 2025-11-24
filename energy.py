"""
Module 
"""
from pathlib import Path
from argparse import ArgumentParser
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from topo_utils import find_contour, isotropic_dilation
from clawpack.clawutil.data import ClawData
from clawpack.visclaw.gridtools import grid_output_2d
from clawpack.pyclaw.solution import Solution
from clawpack.geoclaw import fgout_tools, fgmax_tools
try:
    from config import config
except ImportError:
    config = dict()


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


def avalanche_specific_energy_volume(q, ref_alt, nx, ny, dl, h0=0, g=9.81):
    h, hu, hv, s = q
    hun = hu*nx + hv*ny
    u2 = divide(hu**2 + hv**2, h**2)
    pot = 1/2 * g * (s-ref_alt) * hun
    kin = 1/2 * u2 * hun
    return np.nansum((pot + kin)*dl), np.nansum(hun*dl)


def lake_specific_energy_volume(q, dx, dy, h0=0, g=9.81):
    h, hu, hv, _ = q
    hu2 = divide(hu**2 + hv**2, h)
    pot = 1/2 * g * (h-h0)**2
    kin = 1/2 * hu2
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
    if len(Z.shape) == 2:
        Z = Z[None, ...]
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


def interpolate_avac_fgmax(outdir, x, y, fgno: int=0):
    fg = fgmax_tools.FGmaxGrid()
    fg.outdir = outdir
    fg.read_fgmax_grids_data(fgno, outdir / "fgmax_grids.data")
    fg.read_output()
    X, Y = np.meshgrid(x, y, indexing="ij")
    return RegularGridInterpolator((fg.x, fg.y), fg.h)((X, Y))


def interpolate_avac_fgout(fg, x, y):
    hmax = fg.read_frame(1).h
    for i, t in enumerate(fg.times[1:], start=2):
        hmax = np.maximum(hmax, fg.read_frame(i).h)
    X, Y = np.meshgrid(x, y, indexing="ij")
    return RegularGridInterpolator((fg.x, fg.y), hmax)((X, Y))


def compute_energies_masses(outdir: str="_output",
                            avac_fgmax_fgno: int=0,
                            wave_fgout_fgno: int=1,
                            avac_fgout_fgno: int=1):

    wave_outdir = Path(outdir).expanduser()

    wave_data = Data(wave_outdir)
    avac_outdir = Path(wave_data.setprob.AVAC_outdir.strip("'")).expanduser()

    avac_data = Data(avac_outdir)

    wave_fg = fgout_tools.FGoutGrid(wave_fgout_fgno, wave_outdir)
    wave_fg.read_fgout_grids_data()

    avac_fg = fgout_tools.FGoutGrid(avac_fgout_fgno, avac_outdir)
    avac_fg.read_fgout_grids_data()

    avac_fgmax_file = avac_outdir / f"fgmax{avac_fgmax_fgno:0>4.0f}.txt"
    if avac_fgmax_fgno > 0 and avac_fgmax_file.exists():
        avac_hmax = interpolate_avac_fgmax(avac_outdir, wave_fg.x, wave_fg.y, avac_fgmax_fgno)
    else:
        avac_hmax = interpolate_avac_fgout(avac_fg, wave_fg.x, wave_fg.y)

    wave_fgout_init = wave_fg.read_frame(1)
    avalanche = isotropic_dilation(avac_hmax > 0, 2) & ~isotropic_dilation(wave_fgout_init.h > 0, 2)
    lake_alt = config["lake_alt"]

    avac_fg.extent = avac_fg.x1, avac_fg.x2, avac_fg.y1, avac_fg.y2
    wave_fg.extent = wave_fg.x1, wave_fg.x2, wave_fg.y1, wave_fg.y2
    xc, yc = find_contour(avalanche[:, ::-1], wave_fg.extent, wave_fg.nx, wave_fg.ny).T
    nx, ny, dl = normal_vectors(xc, yc)

    dEa = np.zeros(avac_fg.nout, np.float64)
    dMa = np.zeros(avac_fg.nout, np.float64)
    for i, t in enumerate(avac_fg.times):
        q = avac_fg.read_frame(i+1).q
        qc = np.array((
            RegularGridInterpolator((avac_fg.x, avac_fg.y), q[0])((xc, yc)),
            RegularGridInterpolator((avac_fg.x, avac_fg.y), q[1])((xc, yc)),
            RegularGridInterpolator((avac_fg.x, avac_fg.y), q[2])((xc, yc)),
            RegularGridInterpolator((avac_fg.x, avac_fg.y), q[3])((xc, yc)),
        ))
        dEa[i], dMa[i] = avalanche_specific_energy_volume(qc, lake_alt, nx, ny, dl)

    ta, Ea = trapint(avac_fg.times, avac_data.setprob.rho*dEa)  # TODO: set rho in geoclaw instead of septob in AVAC
    ta, Ma = trapint(avac_fg.times, avac_data.setprob.rho*dMa)

    Ea *= np.sign(Ea.mean())
    Ma *= np.sign(Ma.mean())

    dEw = np.zeros(wave_fg.nout, np.float64)
    dMw = np.zeros(wave_fg.nout, np.float64)
    El = np.zeros(wave_fg.nout, np.float64)
    Ml = np.zeros(wave_fg.nout, np.float64)
    dx = (wave_fg.x2 - wave_fg.x1) / wave_fg.nx
    dy = (wave_fg.y2 - wave_fg.y1) / wave_fg.ny
    xc = np.clip(xc, wave_fg.x.min(), wave_fg.x.max())
    yc = np.clip(yc, wave_fg.y.min(), wave_fg.y.max())
    for i, t in enumerate(wave_fg.times):
        q = wave_fg.read_frame(i+1).q
        qc = np.array((
            RegularGridInterpolator((wave_fg.x, wave_fg.y), q[0])((xc, yc)),
            RegularGridInterpolator((wave_fg.x, wave_fg.y), q[1])((xc, yc)),
            RegularGridInterpolator((wave_fg.x, wave_fg.y), q[2])((xc, yc)),
            RegularGridInterpolator((wave_fg.x, wave_fg.y), q[3])((xc, yc)),
        ))
        dEw[i], dMw[i] = avalanche_specific_energy_volume(qc, lake_alt, nx, ny, dl)
        q[:, avalanche] = 0.
        El[i], Ml[i] = lake_specific_energy_volume(q, dx, dy, wave_fgout_init.h)

    tw, Ew = trapint(wave_fg.times, wave_data.geoclaw.rho*dEw)
    tw, Mw = trapint(wave_fg.times, wave_data.geoclaw.rho*dMw)
    Ew *= np.sign(Ew.mean())
    Mw *= np.sign(Mw.mean())

    El *= wave_data.geoclaw.rho
    Ml *= wave_data.geoclaw.rho
    tl = wave_fg.times

    return (ta, Ea, Ma), (tw, Ew, Mw), (tl, El, Ml)


def plot(ta, Ea, Ma, tw, Ew, Mw, tl, El, Ml):
    fig, axes = plt.subplots(nrows=2, sharex="all")
    axes[0].plot(ta, Ea/1e9, label="AVAC")
    axes[0].plot(tw, Ew/1e9, label="Avalanche")
    axes[0].plot(tl, (El-El[0])/1e9, label="Lake")
    axes[1].plot(ta, Ma/1e3, label="AVAC")
    axes[1].plot(tw, Mw/1e3, label="Avalanche")
    axes[1].plot(tl, (Ml-Ml[0])/1e3, label="Lake")
    for ax in axes:
        ax.legend()
    axes[0].set_ylabel("Energy [GJ]")
    axes[1].set_ylabel("Mass [t]")
    axes[1].set_xlabel("Time [s]")
    return fig, axes


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("-o", "--outdir", type=str, default="_output")
    parser.add_argument("--wave_fgout_fgno", "-w", type=int, default=1)
    parser.add_argument("--avac_fgout_fgno", "-a", type=int, default=config["fgout_fgno"])
    parser.add_argument("--avac_fgmax_fgno", "-m", type=int, default=config["fgmax_fgno"])
    return parser.parse_args()


def main(*args, **kwargs):
    data = compute_energies_masses(*args, **kwargs)
    plot(*[v for d in data for v in d])
    print(f"{data[0][1].max():.2e}, {data[2][1].max():.2e}")
    plt.show()


if __name__ == "__main__":
    main(**parse_args().__dict__)
