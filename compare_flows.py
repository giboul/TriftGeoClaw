from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load
from matplotlib import pyplot as plt
import numpy as np
from clawpack.visclaw import gridtools
from clawpack.pyclaw.solution import Solution
from clawpack.visclaw import gridtools
from matplotlib.animation import FuncAnimation

with open("config.yaml") as file:
    config = safe_load(file)["AVAC"]

parser = ArgumentParser()
parser.add_argument("avid", nargs="?", default="5")
args = parser.parse_args()

x, y = np.loadtxt(Path("TOPM") / "contour_dilated.xy").T
dist = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
dist = np.hstack((0, dist, 2*dist[-1]-dist[-2]))
outdir1 = Path("AVAC") / f"_output{args.avid}"
outdir2 = Path("TSUL") / f"_output{args.avid}"
nfiles = min(len(list(outdir1.glob("fort.q*"))), len(list(outdir2.glob("fort.q*"))))


def extract(i, outdir):
    frame_sol = Solution(i, path=outdir, file_format=config["out_format"])
    q = gridtools.grid_output_2d(
        frame_sol,
        lambda q: q,
        x, y,
        levels = "all",
        return_ma=True
    )
    return q, frame_sol.t

with plt.style.context("bmh"):
   fig, ax = plt.subplots(layout="tight")
   q, t = extract(0, outdir1)
   h, hu, hv, eta = q
   z = eta - h
   zlow = 1.1*z.min()-0.1*z.max()
   eta_steps1 = ax.stairs(eta, dist, baseline=z, fill=True, label="water", color="skyblue")
   z_steps1 = ax.stairs(z, dist, baseline=zlow,
                       label="land", fill=True, color="sienna", lw=1)
   q, t = extract(0, outdir2)
   h, hu, hv, eta = q
   z = eta - h
   zlow = 1.1*z.min()-0.1*z.max()
   eta_steps2 = ax.stairs(eta, dist, baseline=z, label="water", color="r")
   z_steps2 = ax.stairs(z, dist, baseline=zlow, label="land", color="b", lw=1)
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
   (h, hu, hv, eta), t = extract(i, outdir1)
   z = eta-h
   eta_steps1.set_data(eta, dist, z)
   z_steps1.set_data(z, dist, zlow)
   (h, hu, hv, eta), t = extract(i, outdir2)
   z = eta-h
   eta_steps2.set_data(eta, dist, z)
   z_steps2.set_data(z, dist, zlow)
   ax.set_title(title % t)

anim = FuncAnimation(fig, update, nfiles, interval=500)
plt.show()
