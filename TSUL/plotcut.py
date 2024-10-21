from argparse import ArgumentParser
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from params import bounds, out_format
from clawpack.visclaw import colormaps, frametools, geoplot, gridtools
from clawpack.pyclaw import solution


outdir = "_output"
files = list(Path(outdir).glob('fort.q*'))

y = np.linspace(bounds["ymin"], bounds["ymax"])
x = np.linspace(2670500, 2670100)
l = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
l = np.hstack((0, l, 2*l[-1]-l[-2]))


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
    eta_steps = ax.stairs(eta, l, baseline=eta-h, fill=True,
                          label="water", color="skyblue")
    z_steps = ax.stairs(z, l, baseline=1.1*z.min()-0.1*z.max(),
              label="land", fill=True, color="sienna", lw=1)
    title = "Cut @ t=%.2f"
    ax.set_title(title % t)
    # ax.set_aspect("equal")
    ax.legend()
    ax.set_xlabel("Distance [m]")
    ax.set_ylabel("Elevation [MASL]")

def update(frame_num):
    frame_sol = solution.Solution(frame_num, path="_output", file_format=out_format)

    q, t = extract(frame_num)
    h, hu, hv, eta = q

    z = eta-h
    eta_steps.set_data(eta, l, z)
    z_steps.set_data(z, l, z.min()-0.1*(z.max()-z.min()))
    ax.set_title(title % t)
    fig.canvas.draw()

anim = FuncAnimation(fig, update, len(files), interval=500)

parser = ArgumentParser()
parser.add_argument("-f", "--file", action="store_true")
args = parser.parse_args()
if args.file:
    anim.save("stairs.gif", savefig_kwargs=dict(bbox_inches="tight"))
else:
    plt.show()

