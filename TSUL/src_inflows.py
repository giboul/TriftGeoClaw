from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load
from shutil import rmtree
import numpy as np
from clawpack.visclaw import gridtools
from clawpack.pyclaw.solution import Solution
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection

projdir = Path().absolute().parent
with open(projdir / "config.yaml") as file:
    config = safe_load(file)["AVAC"]

parser = ArgumentParser()
parser.add_argument("avid", nargs="?", default="")
parser.add_argument("-p", "--plot", action="store_true")
parser.add_argument("-m", "--movie", action="store_true")
args = parser.parse_args()


inflowdir = projdir / "TSUL" / f"_inflows{args.avid}"
outdir = projdir / "AVAC" / f"_output{args.avid}"

files = list(outdir.glob("fort.q*"))
n = 100
x, y = np.loadtxt(projdir / "TOPM" / "contour_dilated.xy").T
dist = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
dist = np.hstack((0, dist, 2*dist[-1]-dist[-2]))

def extract(i):
    frame_sol = Solution(i, path=outdir, file_format=config["out_format"])
    q = gridtools.grid_output_2d(
        frame_sol,
        lambda q: q,
        x, y,
        levels = "all",
        return_ma=True
    )
    return q, frame_sol.t

def write():
    rmtree(inflowdir, ignore_errors=True)
    nf = len(files)
    Path(inflowdir).mkdir(exist_ok=True)
    times = []
    for ti in range(nf):
        print(f"Saving cut to '{inflowdir}' {ti+1:>{4}}/{nf}...", end="\r")
        q, t = extract(ti)
        times.append(t)
        h, hu, hv, eta = q
        data = np.column_stack((x, y, h, hu, hv))
        path = inflowdir / f"cut{ti:0>{4}}.txt"
        np.savetxt(path, data, comments="")
    np.savetxt(inflowdir / "timing.txt", times)
    print()


def plot(movie):
    with plt.style.context("bmh"):
        fig, (ax1, ax2) = plt.subplots(ncols=2, layout="tight")
        q, t = extract(0)
        h, hu, hv, eta = q
        z = eta - h
        zlow = 1.1*z.min()-0.1*z.max()
        eta_steps = ax1.stairs(eta, dist, baseline=z, fill=True, label="water", color="skyblue")
        z_steps = ax1.stairs(z, dist, baseline=zlow,
                            label="land", fill=True, color="sienna", lw=1)
        title = "Cut @ t=%.2f"
        ax1.set_title(title % t)
        text_z = 0.9*eta.max()+0.1*zlow
        prev_dist = 0
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
        # contour = ax2.scatter(x, y)
    
    def update(i):
        (h, hu, hv, eta), t = extract(i)
        z = eta-h
        eta_steps.set_data(eta, dist, z)
        z_steps.set_data(z, dist, zlow)
        ax1.set_title(title % t)
        lc.set_color(plt.cm.viridis(h))
        lc.set_clim(vmin=h.min(), vmax=h.max())
        # contour.set_color(plt.cm.viridis(h))
        # contour.set_clim(vmin=h.min(), vmax=h.max())

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
