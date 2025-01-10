from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from clawpack.pyclaw.solution import Solution
from clawpack.visclaw.gridtools import grid_output_2d
from maketopo import X, Y, Z, H
from matplotlib.animation import FuncAnimation


outdir = Path("_output")
files = list(outdir.glob("fort.q*"))
x, y = X[X.shape[0]//2, :], Y[X.shape[0]//2, :]


def extract(i):
    frame_sol = Solution(int(i), path=str(outdir), file_format="binary")
    q = grid_output_2d(frame_sol, lambda q: q, x, y, levels = "all", return_ma=True)
    return q, frame_sol.t

def plot():
    # with plt.style.context("dark_background"):
    with plt.style.context("seaborn-v0_8-paper"):
        fig, ax = plt.subplots(layout="tight")
        q, t = extract(1)
        h, hu, hv, eta = q
        z = eta -h
        zlow = z.min()-0.1
        ax.stairs(np.where((np.hstack((h[0], h))[:-1]>1e-10), z+h, np.nan)[:-1], x, baseline=0, label="Surface initiale", fill=False, lw=2)
        land = ax.stairs(z[:-1], x, baseline=zlow, facecolor="sienna", fill=True)
        line = ax.stairs((z+h)[:-1], x, baseline=z[:-1], label=f"Surface à ${t=:.2f}$", fill=True, color="skyblue")
        ax.legend(loc="upper right")
        ax.set_xlabel("$x$ [m]")
        ax.set_ylabel("$z$ [m]")
        ax.set_xlim(x.min(), 2)
        ax.set_ylim(-0.1, 1)
        ax.set_aspect("equal")

    def update(i):
        (h, hu, hv, eta), t = extract(i)
        z = eta - h
        land.set_data(z[:-1], x, -0.1)
        line.set_data((z+h)[:-1], x, z[:-1])
        line.set_label(f"Surface à $h({t=:.2f})$")
        ax.legend(loc="upper right")
        fig.canvas.draw()
    
    return FuncAnimation(fig, update, len(files), interval=500)

anim = plot()
state = dict(v=False)
def pause(e):
    if e.key==" ":
        if state["v"] == True:
            anim.resume()
        else:
            anim.pause()
        state["v"] = not state["v"]
plt.gcf().canvas.mpl_connect("key_press_event", pause)
plt.show()
