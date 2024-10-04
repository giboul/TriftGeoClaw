from pathlib import Path
from params import out_format
from lake_bounds import xmin, xmax, ymin, ymax
import numpy as np
from clawpack.visclaw import gridtools
from clawpack.pyclaw import solution
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


files = list(Path("_output").glob("fort.q*"))
n = 100
x = np.hstack((
    np.linspace(xmin, xmax, n),
    np.full(n, xmax),
    np.linspace(xmax, xmin, n),
    np.full(n, xmin),
))
y = np.hstack((
    np.full(n, ymin),
    np.linspace(ymin, ymax, n),
    np.full(n, ymax),
    np.linspace(ymax, ymin, n)
))

dist = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
dist = np.hstack((0, dist, 2*dist[-1]-dist[-2]))

def extract(i):
    frame_sol = solution.Solution(i, path="_output", file_format=out_format)
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
    q, t = extract(0)
    h, hu, hv, eta = q
    z = eta - h
    eta_steps = ax.stairs(eta, dist, baseline=z, fill=True,
                          label="water", color="skyblue")
    z_steps = ax.stairs(z, dist, baseline=1.1*z.min()-0.1*z.max(),
                        label="land", fill=True, color="sienna", lw=1)
    title = "Cut @ t=%.2f"
    ax.set_title(title % t)
    # ax.set_aspect("equal")
    ax.legend()
    ax.set_xlabel("Distance [m]")
    ax.set_ylabel("Elevation [MASL]")

def update(i):
    (h, hu, hv, eta), t = extract(i)
    z = eta-h
    eta_steps.set_data(eta, dist, z)
    z_steps.set_data(z, dist, 1.1*z.min()-0.1*z.max())
    ax.set_title(title % t)

anim = FuncAnimation(fig, update, len(files), interval=500)
plt.show()

