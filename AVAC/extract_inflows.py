from pathlib import Path
from AddSetrun import out_format
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
boundaries = ("bottom", "right", "top", "left")

dist1 = xmax-xmin
dist2 = dist1 + ymax-ymin
dist3 = dist2 + dist1
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

nf = len(files)
max_digits = len(str(nf))
Path("_cut_output").mkdir(exist_ok=True)
for ti in range(nf):
    print(f"Saving cut {ti+1:>{max_digits}}/{nf}...", end="\r")
    q, t = extract(ti)
    h, hu, hv, eta = q
    for bi, boundary in enumerate(boundaries):
        s = slice(bi*n, (bi+1)*n)
        data = np.vstack((x[s], y[s], h[s], hu[s], hv[s])).T
        np.savetxt(Path("_cut_output") / f"cut.{boundary}_{ti:0>{max_digits}}",
                   data, header=f"{t} := t", comments="")
print()

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
    for d in (dist1, dist2, dist3):
        plt.axline((d, z.min()), slope=float("inf"), ls="-.", c="k")
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

