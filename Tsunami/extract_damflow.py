from pathlib import Path
import numpy as np
from clawpack.visclaw import gridtools
from clawpack.pyclaw import solution
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import params
from maketopo import dam_downstream, dam_upstream


xmin, xmax, ymin, ymax = params.bounds.values()

files = list(Path("_output").glob("fort.q*"))
n = 100
x = np.linspace(xmax, xmin, n, endpoint=True)
y = (dam_downstream(x) + dam_upstream(x))/2
x = x[y <= ymax]
y = y[y <= ymax]
dist = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
dist = np.hstack((0, dist, 2*dist[-1]-dist[-2]))

def extract(i):
    frame_sol = solution.Solution(i, path="_output", file_format=params.out_format)
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
        zlow = 1.1*z.min()-0.1*z.max()
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
        ax.set_ylim(zlow, eta.max())

    def update(i):
        (h, hu, hv, eta), t = extract(i)
        z = eta-h
        eta_steps.set_data(eta, dist, z)
        z_steps.set_data(z, dist, zlow)
        ax.set_title(title % t)

    anim = FuncAnimation(fig, update, len(files), interval=500)
    plt.show()


if __name__ == "__main__":
    plot()
