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

x, y = np.loadtxt(projdir / "TOPM" / "contour_dilated.xy").T
dist = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
dist = centered(np.hstack((0, dist)))
avid = args.avid or ""
outdirAVAC = projdir / "AVAC" / f"_output{avid}"
outdirTSUL = projdir / "TSUL" / f"_output{avid}"
nfiles = min(len(list(outdirAVAC.glob("fort.q*"))), len(list(outdirTSUL.glob("fort.q*"))))


def extract(i, outdir, format=config["out_format"]):
    frame_sol = Solution(i, path=outdir, file_format=format)
    q = gridtools.grid_output_2d(frame_sol, lambda q: q, x, y, levels="all", return_ma=True)
    return q, frame_sol.t

with plt.style.context("bmh"):
   fig, (ax1, ax2) = plt.subplots(layout="tight", ncols=2, sharex=True, sharey=True)
   q, t = extract(0, outdirAVAC)
   h, hu, hv, eta = q
   z = eta - h
   zlow = 1.1*z.min()-0.1*z.max()
   eta_steps1 = ax1.stairs(eta, dist, baseline=z, fill=True, label="water", color="skyblue")
   z_steps1 = ax1.stairs(z, dist, baseline=zlow, label="land", fill=True, color="sienna", lw=1)
   q, t = extract(0, outdirTSUL)
   h, hu, hv, eta = q
   z = eta - h
   zlow = 1.1*z.min()-0.1*z.max()
   eta_steps2 = ax2.stairs(eta, dist, baseline=z, fill=True, label="water", color="skyblue")
   z_steps2 = ax2.stairs(z, dist, baseline=zlow, label="land", fill=True, color="sienna", lw=1)
   title = "Cut @ t=%.2f"
   fig.suptitle(title % t)
   text_z = 0.9*eta.max()+0.1*zlow
   prev_dist = 0
   # ax.set_aspect("equal")
   ax1.legend(loc="lower right")
   ax2.legend(loc="lower right")
   ax1.set_xlabel("Distance [m]")
   ax1.set_ylabel("Elevation [MASL]")
   ax2.set_xlabel("Distance [m]")
   ax2.set_ylabel("Elevation [MASL]")
   ax1.set_xlim(dist[0], dist[-1])
   ax1.set_ylim(zlow, eta.max())
   ax2.yaxis.tick_left()
   ax2.yaxis.set_label_position("right")
   ax1.set_title("AVAC")
   ax2.set_title("TSUL")

def update(i):
   try: i=int(i)
   except ValueError: return None
   (h, hu, hv, eta), t = extract(i, outdirAVAC)
   z = eta-h
   eta_steps1.set_data(eta, dist, z)
   z_steps1.set_data(z, dist, zlow)
   (h, hu, hv, eta), t = extract(i, outdirTSUL)
   z = eta-h
   eta_steps2.set_data(eta, dist, z)
   z_steps2.set_data(z, dist, zlow)
   fig.suptitle(title % t)
   fig.canvas.draw()

fig.canvas.mpl_connect("key_press_event", lambda e: update(e.key))
plt.show()
exit()
anim = FuncAnimation(fig, update, nfiles, interval=500)
plt.show()
