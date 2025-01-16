from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load
import numpy as np
from clawpack.visclaw import gridtools
from clawpack.pyclaw import solution
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    TOPM = config["TOPM"]
    TSUL = config["TSUL"]

bounds = TSUL.get("bounds") or TOPM["bounds"]

parser = ArgumentParser()
parser.add_argument("avid", nargs="?", default="")
parser.add_argument("-m", "--movie", action="store_true")
args = parser.parse_args()
outdir = projdir / "TSUL" / f"_output{args.avid}"
files = list(outdir.glob("fort.q*"))
x, y = np.loadtxt(projdir / "TOPM" / "dam.xy").T
xmin = x[0]-(x[2]-x[0])
x += 20
y -= 20
x = np.hstack((xmin, x))
y = np.hstack((y[0], y))
dist = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
dist = np.hstack((0, dist, 2*dist[-1]-dist[-2]))

def extract(i):
    frame_sol = solution.Solution(i, path=outdir, file_format=TSUL["out_format"])
    q = gridtools.grid_output_2d(
        frame_sol,
        lambda q: q,
        x, y,
        levels = "all",
        return_ma=True
    )
    return q, frame_sol.t

def plot():
    state = dict(zmax=TOPM["lake_alt"])
    with plt.style.context("bmh"):
        fig, ax = plt.subplots(layout="tight")
        q, t = extract(0)
        h, hu, hv, eta = q
        z = eta - h
        zlow = TOPM["dam_alt"]-10
        ax.fill_between(dist, zlow, TOPM["dam_alt"], label="Parapet", facecolor="grey")
        eta_steps = ax.stairs(eta, dist, baseline=z, fill=True, label="Lac", color="skyblue")
        # ax.stairs(eta, dist, baseline=z, fill=True, label="$t=0$", color="royalblue", zorder=eta_steps.get_zorder()+0.1, alpha=0.4)
        z_steps = ax.stairs(z, dist, baseline=zlow,
                            label="Roche", fill=True, color="sienna", lw=1)
        line = ax.axline((x[0], state["zmax"]), slope=0, ls='-.', c="salmon")
        title = "Niveau d'eau au barrage (t=%.2f)"
        ax.set_title(title % t)
        # ax.set_aspect("equal")
        ax.legend(loc="upper right")
        ax.set_xlabel("Distance [m]")
        ax.set_ylabel("Niveau [m.s.m.]")
        ax.set_xlim(dist[0], dist[-1])
        ax.set_ylim(TOPM["dam_alt"]-10, TOPM["dam_alt"]+20)


    def update(i):
        (h, hu, hv, eta), t = extract(i)
        z = eta-h
        state["zmax"] = max(state["zmax"], np.where(h>1e-5, eta, 0).max())
        line.set_xy1((x[0], state["zmax"]))
        eta_steps.set_data(eta, dist, z)
        z_steps.set_data(z, dist, zlow)
        ax.set_title(title % t)

    anim = FuncAnimation(fig, update, len(files), interval=500)
    if args.movie:
        anim.save("cutmovie.gif")
    else:
        plt.show()


if __name__ == "__main__":
    plot()

