#!/usr/bin/env python
from argparse import ArgumentParser
from pathlib import Path
import numpy as np
import pyvista as pv
from clawpack.geoclaw import fgout_tools
from clawpack.clawutil.data import ClawData


def animation(outdir, color_var="dh", fgno=1, cmaps=("qist_earth", "jet"), clim=(-1, 1), file_name=""):

    outdir = Path(outdir)
    clawdata = ClawData()
    clawdata.read(outdir / "claw.data", force=True)
    geodata = ClawData()
    geodata.read(outdir / "geoclaw.data", force=True)
    probdata = ClawData()
    probdata.read(outdir / "setprob.data", force=True)
    fgout_grid = fgout_tools.FGoutGrid(fgno, outdir)
    fgout_grid.read_fgout_grids_data()

    X = fgout_grid.X
    Y = fgout_grid.Y

    fgout_init = fgout_grid.read_frame(1)
    fgout_init.dh = fgout_init.eta - fgout_init.eta

    p = pv.Plotter()

    bathy = pv.StructuredGrid(X, Y, fgout_init.B)
    bathy["z"] = bathy.points[:, 2]
    p.add_mesh(bathy, scalars="z", cmap=cmaps[0], clim=(0, 2500), show_scalar_bar=False)

    surf = pv.StructuredGrid(X, Y, np.where(fgout_init.h>0, fgout_init.eta, np.nan))
    surf[color_var] = get_value(fgout_init, fgout_init, color_var).T.flatten()
    p.add_mesh(surf, scalars=color_var, cmap=cmaps[1], clim=clim, show_scalar_bar=True)

    state = dict(i=0)
    def update(i):
        i += 1
        fgout = fgout_grid.read_frame(i)
        fgout.dh = fgout.eta - fgout_init.eta
        bathy.points[:, 2] = fgout.B.T.flatten()
        surf.points[:, 2] = np.where(fgout.h>0, fgout.eta, np.nan).T.flatten()
        surf[color_var] = get_value(fgout, fgout_init, color_var).T.flatten()

    def update_index(increment=None, value=None):
        if increment is not None:
            state["i"] += increment
        else:
            state["i"] = value
        update(state["i"])
        p.update()

    def prev_frame():
        update_index(increment=-1)
    def next_frame():
        update_index(increment=+1)
    def prevv_frame():
        update_index(increment=-10)
    def nextt_frame():
        update_index(increment=+10)

    if file_name:
        p.open_gif(file_name)
        for i, t in enumerate(fgout_grid.times):
            update(i)
            p.write_frame()
        p.close()
    elif 1:
        p.add_key_event("h", prevv_frame)
        p.add_key_event("j", prev_frame)
        p.add_key_event("k", next_frame)
        p.add_key_event("l", nextt_frame)
        p.show()
    else:
        p.add_timer_event(max_steps=len(fgout_grid.times), duration=500, callback=update)
        p.show()


def get_value(fgout, fgout0, name):
    if name[0] == "d" and hasattr(fgout, name[1:]):
        return getattr(fgout, name[1:]) - getattr(fgout0, name[1:])
    return getattr(fgout, name)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("outdir", type=str, nargs="?", default="_output")
    parser.add_argument("--color_var", "-c", type=str, default="dh")
    parser.add_argument("--fgno", "-n", type=int, default=1)
    parser.add_argument("--cmaps", "-m", type=str, nargs=2, default=("gist_earth", "RdBu"))
    parser.add_argument("--clim", "-l", type=float, nargs=2, default=(-0.5, 0.5))
    parser.add_argument("--file_name", "-f", type=str, default="")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    animation(**args.__dict__)

if __name__ == "__main__":
    main()
