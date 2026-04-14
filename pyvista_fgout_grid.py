#!/usr/bin/env python
from argparse import ArgumentParser
from pathlib import Path
import numpy as np
import pyvista as pv
from matplotlib import colors
from clawpack.geoclaw import fgout_tools


pv.global_theme.allow_empty_mesh = True

def read_asc_dataline(file):
    return file.readline().strip().split()

def read_asc(path, flat=False, dtype=np.float16):
    print(f"{__name__}.read_asc: Reading {path}", end="... ", flush=True)

    path = Path(path)
    if path.suffix == ".asc":
        order = "C"
        datapath = path
        Z = np.loadtxt(path, skiprows=6)
    elif path.suffix == ".bin":
        order = "F"
        datapath = path.with_suffix(".data")
        Z = np.fromfile(path, dtype=dtype)
    else:
        print()
        raise NotImplementedError(f"{path.suffix} not yet supported.")

    with open(datapath) as datafile:
        nx = int(read_asc_dataline(datafile)[0])
        ny = int(read_asc_dataline(datafile)[0])
        xmin = float(read_asc_dataline(datafile)[0])
        ymin = float(read_asc_dataline(datafile)[0])
        resolution = float(read_asc_dataline(datafile)[0])
        _ = float(read_asc_dataline(datafile)[0])
    
    x = xmin + np.arange(nx)*resolution
    y = (ymin + np.arange(ny)*resolution)[::-1]

    print("Done.", flush=True)

    return x, y, Z.reshape(ny, nx, order=order)

def animation(outdir,
              color_var="dh",
              fgno=1,
              cmap="qist_earth",
              cmap_bathy="jet",
              clim=[-1, 1],
              clim_bathy=[0, 2500],
              file_name="",
              init_frame=1,
              animate=False,
              world_png="",
              colorbar_label="",
              hillshade_params=(None, None),
              filter_size=50,
              topo_file="",
              topo_margin=0.0):

    outdir = Path(outdir)
    fgout_grid = fgout_tools.FGoutGrid(fgno, outdir)
    fgout_grid.read_fgout_grids_data()

    X = fgout_grid.X
    Y = fgout_grid.Y

    fgout_init = fgout_grid.read_frame(init_frame)
    fgout_init.dh = fgout_init.eta - fgout_init.eta

    p = pv.Plotter()

    if topo_file:
        x, y, Z = read_asc(topo_file)
        bathy_shape = Z.shape
        Xb, Yb = np.meshgrid(x, y)
        bathy = pv.StructuredGrid(Xb, Yb, Z-topo_margin)
    else:
        bathy_shape = fgout_init.B.shape
        bathy = pv.StructuredGrid(X, Y, fgout_init.B)

    tex = None
    if world_png:
        tex = pv.read_texture(world_png)
        tex.repeat = False
        nx, ny, nc = tex.to_image().dimensions
        
        dx, r1, r2, dy, xmin, ymax = np.loadtxt(Path(world_png).with_suffix(".pgw"))
        xmax = xmin + nx * dx
        ymin = ymax + ny * dy
        o = xmin, ymin, 0
        u = xmax, ymin, 0
        v = xmin, ymax, 0
        bathy.texture_map_to_plane(origin=o, point_u=u, point_v=v, inplace=True)
        scalars = None
    elif hillshade_params[0] is not None:
        ls = colors.LightSource(azdeg=hillshade_params[0], altdeg=hillshade_params[1])
        scalars = "hillshade"
        bathy[scalars] = ls.hillshade(bathy.points[:, 2].reshape(bathy_shape, order="F"), filter_size).flatten(order="F")
        clim_bathy = (0, 1)
        cmap_bathy = "gray"
    else:
        scalars = "z"
        bathy[scalars] = bathy.points[:, 2]

    p.add_mesh(bathy, texture=tex, scalars=scalars, cmap=cmap_bathy, clim=clim_bathy, show_scalar_bar=False)

    surf = pv.StructuredGrid(X, Y, np.where(fgout_init.h>0, fgout_init.eta, np.nan))
    label = colorbar_label or color_var
    surf[label] = get_value(fgout_init, fgout_init, color_var).T.flatten()
    p.add_mesh(surf, scalars=label, cmap=cmap, clim=clim, show_scalar_bar=True)

    state = dict(i=init_frame)
    def update(i):
        fgout = fgout_grid.read_frame(i)
        fgout.dh = fgout.eta - fgout_init.eta
        if not topo_file:
            bathy.points[:, 3] = fgout.B.T.flatten()
        surf.points[:, 2] = np.where(fgout.h>0, fgout.eta, np.nan).T.flatten()
        surf[label] = get_value(fgout, fgout_init, color_var).T.flatten()

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

    def save_movie():
        p.open_gif(file_name)
        s = state["i"]
        for i, t in enumerate(fgout_grid.times[s:], start=s):
            update(i)
            p.write_frame()
    
    def save_image():
        snapshot = "snapshot_"
        i = max([int(path.stem[len(snapshot):]) for path in Path().glob(f"{snapshot}*.svg")] or [0]) + 1
        filename = f"{snapshot}{i}.svg"
        p.save_graphic(
            filename,
            title=f"Frame {state['i']}: t={fgout_grid.times[state['i']-1]}"
        )

    if animate is True:
        p.add_timer_event(max_steps=len(fgout_grid.times), duration=500, callback=update)
        p.show()
    else:
        p.add_key_event("h", prevv_frame)
        p.add_key_event("j", prev_frame)
        p.add_key_event("k", next_frame)
        p.add_key_event("l", nextt_frame)
        p.add_key_event("m", save_movie)
        p.add_key_event("s", save_image)
        p.add_key_event("Left", prev_frame)
        p.add_key_event("Right", next_frame)
        p.show()


def get_value(fgout, fgout0, name):
    if name[0] == "d" and hasattr(fgout, name[1:]):
        return getattr(fgout, name[1:]) - getattr(fgout0, name[1:])
    return getattr(fgout, name)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--outdir", "-o", type=str, default="_output")
    parser.add_argument("--color_var", "-c", type=str, default="dh")
    parser.add_argument("--fgno", "-n", type=int, default=1)
    parser.add_argument("--cmap", "-m", type=str, default="RdBu_r")
    parser.add_argument("--cmap_bathy", type=str, default="gist_earth")
    parser.add_argument("--clim", "-l", type=float, nargs=2, default=[-0.5, 0.5])
    parser.add_argument("--clim_bathy", "-b", type=float, nargs=2, default=[-0.5, 0.5])
    parser.add_argument("--file_name", "-f", type=str, default="fgout.gif")
    parser.add_argument("--init_frame", "-i", type=int, default=1)
    parser.add_argument("--animate", "-a", action="store_true")
    parser.add_argument("--world_png", "-w", type=Path)
    parser.add_argument("--hillshade_params", "-s", type=float, nargs=2, default=[315, 45])
    parser.add_argument("--colorbar_label", "-t", type=str)
    parser.add_argument("--topo_file", type=str, default="")
    args = parser.parse_args()
    return args

def main():
    # args = parse_args()
    # animation(**args.__dict__)
    animation(
        outdir="_output",
        color_var="dh",
        fgno=1,
        cmap="jet",
        cmap_bathy="qist_earth",
        clim=[-1, 1],
        clim_bathy=[0, 2500],
        file_name="movie.gif",
        init_frame=1,
        animate=False,
        world_png="",#"gfx/topo.png",
        colorbar_label="Δh",
        hillshade_params=[315, 45],
        topo_file="",#"topo_trift.asc",
        topo_margin=2  # m
    )

if __name__ == "__main__":
    main()
