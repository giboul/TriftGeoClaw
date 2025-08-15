from argparse import ArgumentParser
from pathlib import Path
import numpy as np
import pyvista as pv
from utils import read_times, read_clawdata, read_fgout_bin


q0 = [0, 0, 0, 0]
def dh(q, q0=q0): return q[0] - q0[0]
def ds(q, q0=q0): return q[3] - q0[3]
def h(q, q0=q0): return q[0]
def hu(q, q0=q0): return q[1]
def hv(q, q0=q0): return q[2]
def b(q, q0=q0): return q[3]-q[0]
def z(q, q0=q0): return q[3]

color_functions = dict(dh=dh, ds=ds, h=h, hu=hu, hv=hv, b=b, z=z)

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("outdir", type=str, nargs="?", default="_output")
    parser.add_argument("--color_by", "-c", type=str, default="dh")
    parser.add_argument("--gridno", "-g", type=int, default=1)
    parser.add_argument("--cmaps", "-m", type=str, nargs=2, default=("gist_earth", "jet"))
    parser.add_argument("--clim", "-l", type=float, nargs=2, default=(-1, 1))
    parser.add_argument("--file_name", "-f", type=str, default="")
    args = parser.parse_args()
    return args

def animation(outdir, color_by="dh", gridno=1, cmaps=("qist_earth", "jet"), clim=(-1, 1), file_name=""):

    outdir = Path(outdir)
    color_func = color_functions[color_by]

    times = read_times(outdir, gridno)

    data = read_clawdata(outdir/"fgout_grids.data", sep="#", skiprows=7)
    x1, y1 = data["x1, y1"]
    x2, y2 = data["x2, y2"]
    nx, ny = data["nx,ny"]
    x = np.linspace(x1, x2, nx)
    y = np.linspace(y1, y2, ny)
    Y, X = np.meshgrid(y, x, copy=False)

    q0 = read_fgout_bin(outdir, nx, ny, frameno=0)

    p = pv.Plotter()

    bathy = pv.StructuredGrid(X, Y, b(q0))
    bathy["z"] = bathy.points[:, 2]
    p.add_mesh(bathy, scalars="z", cmap=cmaps[0], clim=(0, 2500), show_scalar_bar=False)

    surf = pv.StructuredGrid(X, Y, np.where(h(q0)>0, z(q0), np.nan))
    surf[color_by] = color_func(q0, q0).T.flatten()
    p.add_mesh(surf, scalars=color_by, cmap=cmaps[1], clim=clim, show_scalar_bar=True)

    state = dict(i=0)
    def update(i):
        q = read_fgout_bin(outdir, nx, ny, frameno=i)
        bathy.points[:, 2] = b(q).T.flatten()
        surf.points[:, 2] = np.where(h(q)>0, z(q), np.nan).T.flatten()
        surf[color_by] = color_func(q, q0).T.flatten()

    def next_frame():
        state["i"] += 1
        update(state["i"])
        p.update()

    def prev_frame():
        state["i"] -= 1
        update(state["i"])
        p.update()

    if file_name:
        p.open_gif(file_name)
        for i, t in enumerate(times[:20]):
            update(i)
            p.write_frame()
        p.close()
    elif 1:
        p.add_key_event("k", next_frame)
        p.add_key_event("j", prev_frame)
        p.show()
    else:
        p.add_timer_event(max_steps=len(times), duration=500, callback=update)
        p.show()

def main():
    args = parse_args()
    animation(**args.__dict__)

if __name__ == "__main__":
    main()
