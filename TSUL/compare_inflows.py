from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load
from matplotlib import pyplot as plt
import numpy as np
from clawpack.visclaw import gridtools
from clawpack.pyclaw.solution import Solution
from clawpack.visclaw import gridtools
from matplotlib.animation import FuncAnimation

projdir = Path().absolute().parent
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    topoconfig = config["TOPM"]
    config = config["AVAC"]

parser = ArgumentParser()
parser.add_argument("avid", nargs="?", default="")
args = parser.parse_args()

def centered(x):
    x = np.pad(x, (1, 1), mode="edge")
    return (x[1:]+x[:-1])/2

xmin, xmax, ymin, ymax = np.loadtxt(projdir / "TSUL" / "lake_extent.txt")
n = 100
x = np.hstack((
    np.full(n, xmin),
    np.linspace(xmin, xmax, n),
    np.full(n, xmax),
    np.linspace(xmax, xmin, n)
))
y = np.hstack((
    np.linspace(ymax, ymin, n),
    np.full(n, ymin),
    np.linspace(ymin, ymax, n),
    np.full(n, ymax)
))
dist = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
dist = centered(np.hstack((0, dist)))
avid = args.avid or ""
outdirAVAC = projdir / "AVAC" / f"_output{avid}"
outdirTSUL = projdir / "TSUL" / f"_output{avid}"
nfiles = min(len(list(outdirAVAC.glob("fort.q*"))), len(list(outdirTSUL.glob("fort.q*"))))


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
   q, t = extract(0, outdirAVAC)
   h, hu, hv, eta = q
   z = eta - h
   zlow = 1.1*z.min()-0.1*z.max()
   eta_steps1 = ax.stairs(eta, dist, baseline=z, fill=True, label="AVAC: water", color="skyblue")
   z_steps1 = ax.stairs(z, dist, baseline=zlow,
                        label="AVAC: land", fill=True, color="sienna", lw=1)
   q, t = extract(0, outdirTSUL)
   h, hu, hv, eta = q
   z = eta - h
   zlow = 1.1*z.min()-0.1*z.max()
   eta_steps2 = ax.stairs(eta, dist, baseline=z, label="TSUL: water", color="r")
   z_steps2 = ax.stairs(z, dist, baseline=zlow, label="TSUL: land", color="b", lw=1)
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
   try: i=int(i)
   except ValueError: return None
   (h, hu, hv, eta), t = extract(i, outdirAVAC)
   z = eta-h
   eta_steps1.set_data(eta, dist, z)
   z_steps1.set_data(z, dist, zlow)
   zAVAC = z
   (h, hu, hv, eta), t = extract(i, outdirTSUL)
   z = eta-h
   eta_steps2.set_data(eta+zAVAC-z, dist, zAVAC)
   z_steps2.set_data(z, dist, zlow)
   ax.set_title(title % t)
   fig.canvas.draw()

fig.canvas.mpl_connect("key_press_event", lambda e: update(e.key))
plt.show()
exit()
anim = FuncAnimation(fig, update, nfiles, interval=500)
plt.show()
