from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load
import numpy as np
from clawpack.pyclaw.solution import Solution
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection

projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    AVAC = safe_load(file)["AVAC"]

parser = ArgumentParser()
parser.add_argument("avid", nargs="?", default="5")
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

files = list(outdir.glob("fort.q*"))

def extract(i):
    frame_sol = Solution(i, path=outdir, file_format=config["out_format"])
    state = frame_sol.states[0]
    return state.q, state.t

def plot(movie):

    with plt.style.context("bmh"):
        fig, ax = plt.subplots()
        (h, hu, hv, eta), t = extract(0)
        z = eta - h
        bathyim = ax.imshow(z.T[::-1,:])
        him = ax.imshow(np.ma.MaskedArray(h.T[::-1,:], mask=h.T[::-1,:]<=0))
        title = "t=%.2f"
        ax.set_title(title % t)
        # ax.set_aspect("equal")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
    
    def update(i):
        (h, hu, hv, eta), t = extract(i)
        z = eta-h
        bathyim.set_data(z.T[::-1,:])
        h = np.ma.MaskedArray(h.T[::-1,:], mask=h.T[::-1,:]<=0)
        him.set_data(h)
        ax.set_title(title % t)


    files = list(outdir.glob("fort.q*"))
    anim = FuncAnimation(fig, update, len(files), interval=500)
    if movie:
        anim.save("cut_movie.gif")
    else:
        plt.show()


if __name__ == "__main__":
    if args.plot or args.movie or True:
        plot(args.movie)
    else:
        write()
