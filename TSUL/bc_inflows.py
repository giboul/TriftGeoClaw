from argparse import ArgumentParser
from pathlib import Path
from shutil import rmtree
from yaml import safe_load
import numpy as np
from clawpack.visclaw import gridtools
from clawpack.pyclaw.solution import Solution
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


projdir = Path().absolute().parent
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    topoconfig = config["TOPM"]
    config = config["TSUL"]

parser = ArgumentParser()
parser.add_argument("avid", nargs="?", default="")
parser.add_argument("-p", "--plot", action="store_true")
parser.add_argument("-m", "--movie", action="store_true")
args = parser.parse_args()


outTSUL = projdir / "AVAC" / f"_output{args.avid}"
outAVAC = projdir / "AVAC" / f"_output{args.avid}"
inflowdir = projdir / "TSUL" / f"_inflows{args.avid}"

files = list(outTSUL.glob("fort.q*"))
xmin, xmax, ymin, ymax = np.loadtxt("lake_extent.txt")
n = 100
x = np.hstack((
    np.linspace(xmin, xmax, n, endpoint=True),  # South
    np.full(n, xmax),  # East
    np.linspace(xmax, xmin, n, endpoint=True),  # North
    np.full(n, xmin),  # West
))
y = np.hstack((
    np.full(n, ymin),
    np.linspace(ymin, ymax, n),
    np.full(n, ymax),
    np.linspace(ymax, ymin, n)
))
boundaries = ("bottom", "right", "top", "left")

dist1 = x.max()-x.min()
dist2 = dist1 + y.max()-y.min()
dist3 = dist2 + dist1
dist4 = dist3 + dist2 - dist1
dist = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
dist = np.hstack((0, dist, 2*dist[-1]-dist[-2]))

def extract(i, outdir, format):
    frame_sol = Solution(i, path=outdir, file_format=format)
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
        print(f"Saving cut {ti+1:>{4}}/{nf}...", end="\r")
        q, t = extract(ti, outTSUL, config["out_format"])
        times.append(t)
        h, hu, hv, eta = q
        for bi, boundary in enumerate(boundaries):
            s = slice(bi*n, (bi+1)*n)
            data = np.vstack((x[s], y[s], h[s], hu[s], hv[s])).T
            path = inflowdir / f"{boundary}_{ti:0>{4}}.txt"
            np.savetxt(path, data, comments="")
    np.savetxt(inflowdir / "timing.txt", times)
    print()


def plot(movie):
    with plt.style.context("bmh"):
        fig, ax = plt.subplots(layout="tight")
        q, t = extract(0, outTSUL, config["out_format"])
        h, hu, hv, eta = q
        z = eta - h
        zlow = 1.1*z.min()-0.1*z.max()
        eta_steps = ax.stairs(eta, dist, baseline=z, fill=True, label="water", color="skyblue")
        z_steps = ax.stairs(z, dist, baseline=zlow, label="land", fill=True, color="sienna", lw=1)
        title = "Cut @ t=%.2f"
        ax.set_title(title % t)
        text_z = 0.9*eta.max()+0.1*zlow
        prev_dist = 0
        for d, direct in zip((dist1, dist2, dist3, dist4), ("South", "East", "North", "West")):
            plt.text((prev_dist+d)/2, text_z, direct, horizontalalignment='center')
            plt.axline((d, zlow), slope=float("inf"), ls="-.", c="k")
            prev_dist = d
        # ax.set_aspect("equal")
        ax.legend(loc="lower right")
        ax.set_xlabel("Distance [m]")
        ax.set_ylabel("Elevation [MASL]")
        ax.set_xlim(dist[0], dist[-1])
        ax.set_ylim(zlow, eta.max())
    
    def update(i):
        (h, hu, hv, eta), t = extract(i, outTSUL, config["out_format"])
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
