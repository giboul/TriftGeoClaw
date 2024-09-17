from argparse import ArgumentParser
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from params import xmin, xmax, ymin, ymax, out_format
from clawpack.visclaw import colormaps, frametools, geoplot, gridtools
from clawpack.pyclaw import solution


outdir = "_output"
files = list(Path(outdir).glob('fort.q*'))

y = np.linspace(ymin, ymax)
x = np.linspace(2670500, 2670100)
l = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
l = np.hstack((0, l, l[-1]+(l[-1]-l[-2])))

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

fig, ax = plt.subplots()
q, t = extract(0)
h, hu, hv, _1, _2, eta = q
eta_steps = ax.stairs(eta, l, baseline=eta-h, fill=True, label=r"$\eta$", color="skyblue")
ax.stairs(eta-h, l, baseline=eta.min()-0.1*(eta.max()-eta.min()), label=r"$z_b$", fill=True, hatch="//", fc="none", ec="sienna", lw=1)

title = f"Cut @(t=%f, x={x[0]:.0f})"
ax.set_title(title % t)
ax.set_aspect("equal")
ax.legend()

def update(frame_num):
    frame_sol = solution.Solution(frame_num, path="_output", file_format=out_format)

    q, t = extract(frame_num)
    h, hu, hv, _1, _2, eta = q

    eta_steps.set_data(eta, l, eta-h)
    ax.set_title(title % t)

anim = FuncAnimation(fig, update, len(files), interval=500)
plt.show()

parser = ArgumentParser()
parser.add_argument("-f", "--file", action="store_true")
args = parser.parse_args()
if args.file:
    html = anim.to_jshtml()
    with open("stairs.html", "w") as file:
        file.write(html)

