from argparse import ArgumentParser
from pathlib import Path
from utils import config, read_world_image, water_color
import topo_utils
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from skimage.morphology import flood, isotropic_dilation
from skimage.measure import find_contours
import matplotlib as mpl
from matplotlib import pyplot as plt
import pyvista as pv


pv.global_theme.allow_empty_mesh = True

mpl.colormaps.register(mpl.colors.LinearSegmentedColormap("half_blue", {
    'red': ((0.0, 1.0, 1.0),
            (1.0, water_color[0], water_color[0]),),
    'green': ((0.0, 0.0, 0.0), 
              (1.0, water_color[1], water_color[1]),),
    'blue': ((0.0, 0.0, 0.0),
             (1.0, water_color[1], water_color[1]),),
    'alpha': ((0.0, 0.0, 0.0),
              (0.99, 0.0, water_color[2]),
              (1.0, water_color[2], water_color[2]),),
}))


def flood_mask(topo, seed, max_level=0):
    mask = topo < max_level
    initial_value = mask[*seed]
    mask[*seed] = True
    flooded = flood(mask, seed)
    flooded[*seed] = initial_value
    return flooded


def iflood(bathy_path, lake_alt=0, dilation_radius=0, world_png=""):

    x, y, Z = topo_utils.read_asc(Path(bathy_path).expanduser())
    Zmin = Z.min()

    pl = pv.Plotter()

    X, Y = np.meshgrid(x, y, copy=False)
    s = (slice(None, None, 1),)*2
    bathy = pv.StructuredGrid(X[s], Y[s], Z[s])
    zmin = bathy.points[:, 2].min()
    water = pv.StructuredGrid(X[s], Y[s], np.full_like(Z[s], zmin))
    water["wet"] = np.zeros(Z[s].size)
    if world_png != "":
        I, extent = read_world_image(world_png)
        world_png = pv.read_texture(Path(world_png).expanduser())
        bathy = bathy.clip_box(extent+(0, 9000), invert=False)
        o = extent[0], extent[2], 0.
        u = extent[1], extent[2], 0.
        v = extent[0], extent[3], 0.
        bathy.texture_map_to_plane(o, u, v, inplace=True)

    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()
    ny, nx = Z[s].shape
    res_x = (xmax - xmin) / nx
    res_y = (ymax - ymin) / ny

    def on_click(point):
        x, y, z = point
        ix = int((point[0]-xmin)/res_x)
        iy = int((ymax-point[1])/res_y)
        ix = max(0, min(Z[s].shape[1], ix))
        iy = max(0, min(Z[s].shape[0], iy))
        water["wet"] = flood_mask(Z[s], (iy, ix), lake_alt).flatten(order="F")
        water["wet"] = isotropic_dilation(water["wet"], radius=dilation_radius)
        water.points[:, 2] = np.where(water["wet"], lake_alt, Zmin)

    pl = pv.Plotter()
    pl.add_mesh(bathy, texture=world_png or None, show_scalar_bar=False)
    pl.add_mesh(water, cmap="half_blue", show_scalar_bar=False)
    pl.enable_surface_point_picking(callback=on_click, show_point=True)
    pl.show()

    print(f"Saving qinit.xyz", end=" ... ", flush=True)
    np.savetxt("qinit.xyz", water.points.T.reshape(3, nx, ny).reshape(3, -1, order="F").T)
    print("Done")

    xc, yc = find_contours(water["wet"].reshape(Z[s].shape, order="F").T, 0.5)[0].T
    xc = xmin + xc * res_x
    yc = ymax - yc * res_y
    np.savetxt("contour.xy", np.column_stack((xc, yc)))


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("bathy_path", default=config.get("bathymetry", "bathymetry.asc"), type=str, nargs="?")
    parser.add_argument("lake_alt", default=config.get("lake_alt", 0.), type=int, nargs="?")
    parser.add_argument("dilation_radius", default=config.get("flood_dilation", 0), type=int, nargs="?")
    parser.add_argument("world_png", default=config.get("world_png", ""), type=str, nargs="?")
    return parser.parse_args()

def main():
    iflood(**parse_args().__dict__)

if __name__ == "__main__":
    main()
