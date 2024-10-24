from argparse import ArgumentParser
from pathlib import Path
from shutil import rmtree
from AddSetrun import out_format, resolution
from lake_bounds import xmin, xmax, ymin, ymax
import numpy as np
from clawpack.visclaw import gridtools
from clawpack.pyclaw import solution
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


parser = ArgumentParser()
parser.add_argument("avid", nargs="?", default="")
parser.add_argument("-p", "--plot", action="store_true")
parser.add_argument("-m", "--movie", action="store_true")
args = parser.parse_args()


outcutdir = Path(f"_cut_output{args.avid}")
outdir = Path(f"_output{args.avid}")

files = list(outdir.glob("fort.q*"))
n = 100
x, y = np.loadtxt("contour_dilated.xy").T

dist = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
dist = np.hstack((0, dist, 2*dist[-1]-dist[-2]))

def extract(i):
    frame_sol = solution.Solution(i, path=outdir, file_format=out_format)
    q = gridtools.grid_output_2d(
        frame_sol,
        lambda q: q,
        x, y,
        levels = "all",
        return_ma=True
    )
    return q, frame_sol.t

def write():
    rmtree(outcutdir, ignore_errors=True)
    nf = len(files)
    Path(outcutdir).mkdir(exist_ok=True)
    times = []
    for ti in range(nf):
        print(f"Saving cut to '{outcutdir}' {ti+1:>{4}}/{nf}...", end="\r")
        q, t = extract(ti)
        times.append(t)
        h, hu, hv, eta = q
        data = np.vstack((x, y, h, hu, hv)).T
        path = outcutdir / f"cut{ti:0>{4}}.txt"
        np.savetxt(path, data, comments="")
    np.savetxt(outcutdir / "timing.txt", times)
    print()


def plot(movie):
    with plt.style.context("bmh"):
        fig, ax = plt.subplots(layout="tight")
        q, t = extract(0)
        h, hu, hv, eta = q
        z = eta - h
        zlow = 1.1*z.min()-0.1*z.max()
        eta_steps = ax.stairs(eta, dist, baseline=z, fill=True, label="water", color="skyblue")
        z_steps = ax.stairs(z, dist, baseline=zlow,
                            label="land", fill=True, color="sienna", lw=1)
        title = "Cut @ t=%.2f"
        ax.set_title(title % t)
        text_z = 0.9*eta.max()+0.1*zlow
        prev_dist = 0
        # ax.set_aspect("equal")
        ax.legend(loc="lower right")
        ax.set_xlabel("Distance [m]")
        ax.set_ylabel("Elevation [MASL]")
        ax.set_xlim(dist[0], dist[-1])
        ax.set_ylim(zlow, eta.max())
    
    def update(i):
        (h, hu, hv, eta), t = extract(i)
        z = eta-h
        eta_steps.set_data(eta, dist, z)
        z_steps.set_data(z, dist, zlow)
        ax.set_title(title % t)

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
