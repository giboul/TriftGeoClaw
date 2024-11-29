from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load
import numpy as np
from clawpack.visclaw import gridtools
from clawpack.pyclaw.solution import Solution
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection

projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    AVAC = safe_load(file)["AVAC"]

parser = ArgumentParser()
parser.add_argument("avid", nargs="?", default="")
parser.add_argument("-p", "--plot", action="store_true")
parser.add_argument("-m", "--movie", action="store_true")
args = parser.parse_args()

outdir = projdir / "TSUL" / f"_output{args.avid}"
AVACoutdir = projdir / "AVAC" / f"_output{args.avid}"

def read_contour(path):
    x, y = np.loadtxt(path).T
    dist = np.cumsum(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2))
    dist = np.hstack((0, dist, 2 * dist[-1] - dist[-2]))
    return np.column_stack((x, y))

contour2 = read_contour(projdir/"TOPM"/"contour2.xy")
contour3 = read_contour(projdir/"TOPM"/"contour3.xy")

def extract(i, x, y):
    frame_sol = Solution(i, path=AVACoutdir, file_format=AVAC["out_format"])
    q = gridtools.grid_output_2d(
        frame_sol, lambda q: q, x, y, levels="all", return_ma=True
    )
    return q, frame_sol.t


def write():

    Path(outdir).mkdir(exist_ok=True)
    for f in outdir.glob("cut*.txt"):
        f.unlink()

    files = list(AVACoutdir.glob("fort.q*"))
    nf = len(files)
    times = []
    for ti in range(nf):
        print(f"Saving cuts to '{outdir}' {ti+1:>{4}}/{nf}...", end="\r")
        for i, contour in enumerate((contour2, contour3), start=2):
            q, t = extract(ti, *contour.T)
            h, hu, hv, eta = q
            data = np.column_stack((*contour.T, h, hu, hv))
            path = outdir / f"cut{i}_{ti:0>{4}}.npy"
            data[~np.isfinite(data)] = 0
            np.savetxt(path, data, comments="")
        times.append(t)
    np.savetxt(outdir / "times.txt", times)
    print()


def plot(movie):

    with plt.style.context("bmh"):
        fig, (ax1, ax2) = plt.subplots(ncols=2, layout="tight")
        x, y = contour2.T
        dist = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
        dist = np.hstack((dist[0], dist, dist[-1]))
        q, t = extract(0, x, y)
        h, hu, hv, eta = q
        z = eta - h
        zlow = 1.1 * z.min() - 0.1 * z.max()
        eta_steps = ax1.stairs(
            eta, dist, baseline=z, fill=True, label="water", color="skyblue"
        )
        z_steps = ax1.stairs(
            z, dist, baseline=zlow, label="land", fill=True, color="sienna", lw=1
        )
        title = "Cut @ t=%.2f"
        ax1.set_title(title % t)
        # ax1.set_aspect("equal")
        ax1.legend(loc="lower right")
        ax1.set_xlabel("Distance [m]")
        ax1.set_ylabel("Elevation [MASL]")
        ax1.set_xlim(dist[0], dist[-1])
        ax1.set_ylim(zlow, eta.max())

        lines = np.c_[x[:-1], y[:-1], x[1:], y[1:]]
        lc = LineCollection(lines.reshape(-1, 2, 2), array=h, linewidths=5)
        ax2.add_collection(lc)
        ax2.set_xlim(x.min(), x.max())
        ax2.set_ylim(y.min(), y.max())

    def update(i):
        x, y = contour2.T
        (h, hu, hv, eta), t = extract(i, x, y)
        z = eta - h
        eta_steps.set_data(eta, dist, z)
        z_steps.set_data(z, dist, zlow)
        ax1.set_title(title % t)
        lc.set_color(plt.cm.viridis(h))
        lc.set_clim(vmin=h.min(), vmax=h.max())


    files = list(outdir.glob("fort.q*"))
    anim = FuncAnimation(fig, update, len(files), interval=500)
    if movie:
        anim.save("cut_movie.gif")
    else:
        plt.show()


if __name__ == "__main__":
    if args.plot or args.movie:
        plot(args.movie)
    else:
        write()