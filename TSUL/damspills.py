from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load
import numpy as np
from clawpack.visclaw import gridtools
from clawpack.pyclaw import solution
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

projdir = Path().absolute().parent
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    topoconfig = config["TOPM"]
    config = config["TSUL"]

bounds = config.get("bounds") or topoconfig["bounds"]

parser = ArgumentParser()
parser.add_argument("avid", nargs="?", default="")
parser.add_argument("-m", "--movie", action="store_true")
args = parser.parse_args()
print(f"{args.avid = }")
outdir = projdir / "TSUL" / f"_output{args.avid}"
print(f"{outdir = }")
files = list(outdir.glob("fort.q*"))
print(f"{len(files) = }")
x, y = np.loadtxt(projdir / "TOPM" / "dam.xy").T
dist = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
dist = np.hstack((0, dist, 2*dist[-1]-dist[-2]))

def extract(i):
    frame_sol = solution.Solution(i, path=outdir, file_format=config["out_format"])
    q = gridtools.grid_output_2d(
        frame_sol,
        lambda q: q,
        x, y,
        levels = "all",
        return_ma=True
    )
    return q, frame_sol.t

def plot():
    with plt.style.context("bmh"):
        fig, ax = plt.subplots(layout="tight")
        q, t = extract(0)
        h, hu, hv, eta = q
        z = eta - h
        zlow = topoconfig["dam_alt"]-10
        eta_steps = ax.stairs(eta, dist, baseline=z, fill=True, label="water", color="skyblue")
        z_steps = ax.stairs(z, dist, baseline=zlow,
                            label="land", fill=True, color="sienna", lw=1)
        title = "Cut @ t=%.2f"
        ax.set_title(title % t)
        # ax.set_aspect("equal")
        ax.legend(loc="upper right")
        ax.set_xlabel("Distance [m]")
        ax.set_ylabel("Elevation [MASL]")
        ax.set_xlim(dist[0], dist[-1])
        ax.set_ylim(topoconfig["dam_alt"]-10, topoconfig["dam_alt"]+20)

    def update(i):
        (h, hu, hv, eta), t = extract(i)
        (h, hu, hv, eta), t = extract(i)
        z = eta-h
        eta_steps.set_data(eta, dist, z)
        z_steps.set_data(z, dist, zlow)
        ax.set_title(title % t)

    anim = FuncAnimation(fig, update, len(files), interval=500)
    if args.movie:
        anim.save("cutmovie.gif")
    else:
        plt.show()


if __name__ == "__main__":
    plot()

