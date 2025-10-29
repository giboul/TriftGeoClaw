#!/usr/bin/env python
from argparse import ArgumentParser
from pathlib import Path
import numpy as np
import pyvista as pv
from clawpack.geoclaw import fgmax_tools
from clawpack.clawutil.data import ClawData


pv.global_theme.allow_empty_mesh = True


def view3d(outdir, color_var="dh", fgno=1, cmaps=("qist_earth", "jet"), clim=(-1, 1), file_name=""):

    outdir = Path(outdir)
    clawdata = ClawData()
    clawdata.read(outdir / "claw.data", force=True)
    geodata = ClawData()
    geodata.read(outdir / "geoclaw.data", force=True)
    probdata = ClawData()
    probdata.read(outdir / "setprob.data", force=True)

    fg = fgmax_tools.FGmaxGrid()
    fg.outdir = outdir
    fg.read_fgmax_grids_data(fgno, outdir / "fgmax_grids.data")
    fg.read_output()

    X = fg.X
    Y = fg.Y

    p = pv.Plotter()

    bathy = pv.StructuredGrid(X, Y, fg.B)
    bathy["z"] = bathy.points[:, 2]
    p.add_mesh(bathy, scalars="z", cmap=cmaps[0], clim=(0, 2500), show_scalar_bar=False)

    surf = pv.StructuredGrid(X, Y, np.where(fg.h>0, fg.B+fg.h, np.nan))
    surf[color_var] = getattr(fg, color_var).T.flatten()
    p.add_mesh(surf, scalars=color_var, cmap=cmaps[1], clim=clim, show_scalar_bar=True)

    p.show()


def get_value(fgmax, fgmax0, name):
    if name[0] == "d" and hasattr(fgmax, name[1:]):
        return getattr(fgmax, name[1:]) - getattr(fgmax0, name[1:])
    return getattr(fgmax, name)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("outdir", type=str, nargs="?", default="_output")
    parser.add_argument("--color_var", "-c", type=str, default="h")
    parser.add_argument("--fgno", "-n", type=int, default=1)
    parser.add_argument("--cmaps", "-m", type=str, nargs=2, default=("gist_earth", "RdBu"))
    parser.add_argument("--clim", "-l", type=float, nargs=2, default=(-0.5, 0.5))
    parser.add_argument("--file_name", "-f", type=str, default="")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    view3d(**args.__dict__)

if __name__ == "__main__":
    main()
