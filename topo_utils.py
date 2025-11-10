#!/usr/bin/env python
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from tifffile import TiffFile
from skimage.morphology import flood
from skimage.measure import find_contours
import geopandas as gpd
import shapely
import rasterio
from clawpack.geoclaw import topotools


def read_raster(path: Path, bbox=None, cellsize=None):
    with rasterio.open(path) as file:
        if bbox is not None:
            w, s, e, n = bbox
            bbox = rasterio.windows.from_bounds(*bbox, transform=file.transform)
        else:
            w, s, e, n = file.bounds
        if cellsize is None:
            out_shape = None
            nx = file.width
            ny = file.height
        else:
            nx = int((e-w)/cellsize)
            ny = int((n-s)/cellsize)
            out_shape = file.count, ny, nx
        Z = file.read(window=bbox, out_shape=out_shape)[0]
    x = np.linspace(w, e, num=nx, endpoint=True)
    y = np.linspace(n, s, num=ny, endpoint=True)
    return x, y, Z


def write_raster(path: Path, Z, bbox, i=1):
    height, width = Z.shape
    w, s, e, n = bbox
    csx = (e - w)/width
    csy = (n - s)/height
    transform = rasterio.transform.Affine(csx, 0, w, 0, -csy, n)
    with rasterio.open(path, "w",
                       height=height,
                       width=width,
                       dtype=Z.dtype,
                       count=1,
                       driver="GTiff",
                       transform=transform) as file:
        file.write(Z, indexes=i)


def read_asc_dataline(file):
    return file.readline().strip().split()


def read_asc(path, dtype=np.float16):
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


def write_geoclaw_topo(path, Z, x, y, topo_type=3):
    topo = topotools.Topography()
    topo.x = x
    topo.y = y
    topo.Z = Z
    topo.write(path, topo_type=topo_type)


def _write_asc(path, Z, xmin, ymin, resolution, nodatavalue=999999, fmt="%.16e", dtype=np.float16):
    path = Path(path)
    print(f"{__name__}.write_asc: Writing to {path}", end="... ", flush=True)
    asc_header = f"""
    {Z.shape[1]} ncols
    {Z.shape[0]} nrows
    {xmin} xllcenter
    {ymin} yllcenter
    {resolution} cellsize
    {nodatavalue} nodata_value
    """.strip()
    if path.suffix == ".asc":
        np.savetxt(path, Z.flatten(), header=asc_header, comments="")
    elif path.suffix == ".bin":
        with open(path.with_suffix(".data"), "w") as file:
            file.write(asc_header)
        Z.reshape(Z.shape, order="F").astype(dtype).T.tofile(path)
    else:
        print()
        raise NotImplementedError(f"{path.suffix} not yet implemented.")
    print("Done.", flush=True)


def read_tiff(path):
    print(f"{__name__}.read_tiff: Reading {path}", end="... ", flush=True)
    with TiffFile(path) as tif:
        Z = tif.asarray().astype(np.float16)
        nx = tif.pages[0].tags["ImageWidth"].value
        ny = tif.pages[0].tags["ImageLength"].value
        x_res, y_res, z_res = tif.pages[0].tags["ModelPixelScaleTag"].value
        xmin = tif.pages[0].tags["ModelTiepointTag"].value[3]
        ymax = tif.pages[0].tags["ModelTiepointTag"].value[4]
        _ = tif.pages[0].tags["ModelPixelScaleTag"].value[0]
    x = xmin + (np.arange(nx) + 0.5)*x_res
    y = ymax - (np.arange(ny) + 0.5)*y_res
    print("Done.")
    return x, y, Z

def read_poly(path):
    geom = gpd.read_file(path).geometry
    if len(geom) == 0:
        raise ValueError(f"No polygon was found in {path}")
    if isinstance(geom[0], shapely.MultiPolygon):
        return [np.array(g.exterior.coords) for poly in geom for g in poly.geoms]
    return [np.array(poly.exterior.coords) for poly in geom]

def flood_mask(topo, seed, max_level=0):
    mask = topo < max_level
    initial_value = mask[*seed]
    mask[*seed] = True
    flooded = flood(mask, seed)
    flooded[*seed] = initial_value
    return flooded


def pick_seed(z_im, x, y, lake_alt=0):

    res = abs(x[1]-x[0])
    extent = x.min(), x.max(), y.min(), y.max()
    fig, ax = plt.subplots(layout="tight")
    imb = ax.imshow(z_im, extent=extent)
    imd = ax.imshow(np.ma.MaskedArray([[1]], mask=True), extent=extent, cmap="Reds")
    imf = ax.imshow(np.ma.MaskedArray([[1]], mask=True), extent=extent, cmap="Blues")
    title = "Max altitude: %.1f   Dilation radius: %i   Status: %s"
    ax.legend(
        [Line2D([0], [0], color=imb.get_cmap()(0.5), lw=4),
         Line2D([0], [0], color=imf.get_cmap()(0.), lw=4),
         Line2D([0], [0], color=imd.get_cmap()(1.), lw=4)],
        ("Bathymetry", "Flooded region", "Dilated flood")
    )

    data = dict(alt=lake_alt, x=0, y=0, r=0, status="waiting")
    dilation_keys = dict(right=+1, left=-1)
    altitude_keys = dict(up=+1, down=-1)

    def on_submit(expression):
        try:
            data["alt"] = float(expression)
        except Exception as e:
            print(e)

    flooded = np.full(z_im.shape, False, dtype=bool)
    dilated = flooded.copy()

    def redraw(ignore_pause=False):
        fig.canvas.manager.set_window_title("Flooding...")
        flooded[:,:] = flood_mask(z_im, (data["y"], data["x"]), float(data.get("alt") or 0))
        dilated[:,:] = isotropic_dilation(flooded, data["r"])
        imf.set_data(np.where(flooded, z_im, np.nan))
        imd.set_data(np.where(dilated, z_im, np.nan))
        if flooded.sum() != 0:
            imf.set(clim=(z_im[flooded].min(), z_im[flooded].max()))
            imd.set(clim=(z_im[dilated].min(), z_im[dilated].max()))
        fig.canvas.draw()
        fig.canvas.manager.set_window_title(f"Status: {data['status']}")

    def key_events(event, draw=True):
        if event.key in dilation_keys:
            data["r"] += dilation_keys[event.key]
        if event.key in altitude_keys:
            data["alt"] += altitude_keys[event.key]
        elif event.key == "enter":
            redraw()
        elif event.key == " ":
            if data["status"] == "pause":
                data["status"] = "waiting"
                redraw()
            else:
                data["status"] = "pause"
        fig.canvas.manager.set_window_title(title % (data["alt"], data["r"], data["status"]))

    fig.canvas.mpl_connect("key_press_event", key_events)

    def flood_pick(event):
        if event.dblclick:
            data["x"] = np.abs(x-event.xdata+res/2).argmin()
            data["y"] = np.abs(y-event.ydata+res/2).argmin()
            redraw(ignore_pause=True)
    fig.canvas.mpl_connect("button_press_event", flood_pick)
    
    plt.show()
    return flooded, dilated


def uniform_grid_interp(x, y, Z, xZ=None, yZ=None):
    if xZ is not None:
        x = (x-xZ[0])/(xZ[-1]-xZ[0])*(Z.shape[1]-1)
    if yZ is not None:
        y = (y-yZ[0])/(yZ[-1]-yZ[0])*(Z.shape[0]-1)
    x = np.clip(0., Z.shape[1], x)
    y = np.clip(0., Z.shape[0], y)
    w = np.clip(x.astype(int), 0, Z.shape[1]-2)
    s = np.clip(y.astype(int), 0, Z.shape[0]-2)
    e = w + 1
    n = s + 1
    return (
        + (x-w)*((y-s)*Z[:, e][n, :].T).T
        + (x-w)*((n-y)*Z[:, e][s, :].T).T
        + (e-x)*((y-s)*Z[:, w][n, :].T).T
        + (e-x)*((n-y)*Z[:, w][s, :].T).T
    )


def read_world_image(path):
    path = Path(path)
    im = plt.imread(path)
    ny, nx, _ = im.shape
    world_text = path.with_suffix(".pgw").read_text().strip().split("\n")
    dx, _, _, dy, xul, yul = [float(f) for f in world_text]
    return im, (xul, xul+nx*dx, yul+ny*dy, yul)


def find_contour(mask, extent, nx, ny):
    xmin, xmax, ymin, ymax = extent
    x, y = find_contours(mask, 0.5)[0].T
    x = xmin + x/nx * (xmax - xmin)
    y = ymax - y/ny * (ymax - ymin)
    return np.column_stack((x, y))


def expand_bounds(x1, x2, y1, y2, rel_margin=1/50, abs_margin=0):
    dx = (x2 - x1) * rel_margin + abs_margin
    dy = (y2 - y1) * rel_margin + abs_margin
    xmin = x1 - dx
    xmax = x2 + dx
    ymin = y1 - dy
    ymax = y2 + dy
    return xmin, xmax, ymin, ymax


def gkern(n, sig):
    """creates gaussian kernel with side length `n` and a sigma of `sig`"""
    ax = np.linspace(-(n - 1) / 2., (n - 1) / 2., n)
    gauss = np.exp(-0.5 * np.square(ax) / np.square(sig))
    kernel = np.outer(gauss, gauss)
    return kernel / np.sum(kernel)


def conv2d(Z, k):
    k = np.array(k)
    if k.size == 0:
        return Z
    k = k / np.sum(k)
    s = k.shape + tuple(np.subtract(Z.shape, k.shape) + 1)
    strd = np.lib.stride_tricks.as_strided
    subM = strd(Z, shape = s, strides = Z.strides * 2)
    return np.einsum('ij,ijkl->kl', k, subM)


def pad_to_shape(Z, shape):
    y_pad, x_pad = np.subtract(shape, Z.shape)
    return np.pad(Z,(
        (y_pad//2, y_pad//2 + y_pad%2), 
        (x_pad//2, x_pad//2 + x_pad%2)
    ), mode = 'edge')


def apply_gkernel(Z, size=1, sig=1.):
    return pad_to_shape(conv2d(Z, gkern(size, sig)), Z.shape)


def isodil(M):
    M[1:-1, 1:-1] = (
        M[:-2,  :-2] | M[1:-1,  :-2] | M[2:,  :-2] |
        M[:-2, 1:-1] | M[1:-1, 1:-1] | M[2:, 1:-1] |
        M[:-2,   2:] | M[1:-1,   2:] | M[2:,   2:]
    )
    return M


def isotropic_dilation(M, r=1):
    for _ in range(r):
        M = isodil(M)
    return M
