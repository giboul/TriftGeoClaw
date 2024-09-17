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
l = np.hstack((0, l))

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
eta_poly = ax.fill_between(l, eta-h, eta, label=r"$\eta$", ec="skyblue")
zb_poly  = ax.fill_between(l, eta-h, eta.min(), label=r"$z_b$", hatch="//", fc="none", ec="sienna")
title = f"Cut @(t=%f, x={x[0]:.0f})"
ax.set_title(title % t)
ax.set_aspect("equal")
ax.legend()

def vertices_between(x, y1, y2):
    new_x = np.hstack((x, x[::-1]))
    new_y = np.hstack((y1, y2[::-1]))
    return np.vstack((new_x, new_y)).T

def update(frame_num):
    frame_sol = solution.Solution(frame_num, path="_output", file_format=out_format)

    q, t = extract(frame_num)
    h, hu, hv, _1, _2, eta = q

    eta_poly.set_verts([vertices_between(l, eta, eta-h)], closed=False)
    ax.set_title(title % t)

anim = FuncAnimation(fig, update, len(files), interval=500)
plt.show()

