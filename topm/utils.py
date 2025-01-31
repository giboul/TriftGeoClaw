#!/usr/bin/env python
from json import load
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import TextBox
from matplotlib.lines import Line2D
from tifffile import TiffFile
from skimage.morphology import flood
from skimage.measure import find_contours


def read_tiff(path):
    with TiffFile(path) as tif:
        Z = tif.asarray().astype(np.float16)
        nx = tif.pages[0].tags["ImageWidth"].value
        ny = tif.pages[0].tags["ImageLength"].value
        x_res, y_res, z_res = tif.pages[0].tags["ModelPixelScaleTag"].value
        xmin = tif.pages[0].tags["ModelTiepointTag"].value[3]
        ymax = tif.pages[0].tags["ModelTiepointTag"].value[4]
        tifres = tif.pages[0].tags["ModelPixelScaleTag"].value[0]
    x = xmin + (np.arange(nx) + 0.5)*x_res
    y = ymax - (np.arange(ny) + 0.5)*y_res
    return x, y, Z

def fill_lake(topo, seed, max_level=0):
    mask = topo < max_level
    initial_value = mask[*seed]
    mask[*seed] = True
    flooded = flood(mask, seed)
    flooded[*seed] = initial_value
    # topo[flooded] = max_level
    return flooded

def read_geojson(path):
    with open(path, "r") as file:
        data = load(file)

    coords = []
    for e in data["features"]:
        ix = e['properties']['id']
        points = e['geometry']['coordinates'][0][0]
        for x, y in points:
            coords.append([ix, x, y])

    avacs = np.array(coords).T
    return avacs

def pick_seed(z_im, x, y, res=0, lake_alt=0):

    extent = x.min(), x.max(), y.min(), y.max()
    fig, (axbox, ax) = plt.subplots(nrows=2, gridspec_kw=dict(height_ratios=(1, 5)), layout="tight")
    text_box = TextBox(axbox, "Lake elevation: ", textalignment="center")
    text_box.set_val(lake_alt)
    imb = ax.imshow(z_im, extent=extent)
    imd = ax.imshow(np.ma.MaskedArray([[1]], mask=True), extent=extent, cmap="Reds")
    imf = ax.imshow(np.ma.MaskedArray([[1]], mask=True), extent=extent, cmap="Blues")
    title = "Max altitude: %s   Dilation radius: %i"
    ax.legend(
        [Line2D([0], [0], color=imb.get_cmap()(0.5), lw=4),
         Line2D([0], [0], color=imf.get_cmap()(0.), lw=4),
         Line2D([0], [0], color=imd.get_cmap()(1.), lw=4)],
        ("Bathymetry", "Flooded region", "Dilated flood")
    )

    data = dict(alt=str(lake_alt), x=0, y=0, r=0, status="waiting")
    keys = dict(up=1, right=1, down=-1, left=-1)
    ax.set_title(title % (data["alt"], data["r"]))

    def on_submit(expression):
        try:
            data["alt"] = float(expression)
        except Exception as e:
            print(e)
    text_box.on_submit(on_submit)

    def redraw(ignore_pause=False):
        if data["status"] == "pause" and ignore_pause is False:
            return None
        fig.canvas.manager.set_window_title("Flooding...")
        flooded = fill_lake(z_im.copy(), (data["y"], data["x"]), float(data.get("alt") or 0))
        dilated = isotropic_dilation(flooded, data["r"])
        flooded = np.ma.MaskedArray(z_im, ~flooded)
        dilated = np.ma.MaskedArray(z_im, ~dilated)
        imf.set_data(flooded)
        imd.set_data(dilated)
        imf.set(clim=(flooded.min(), flooded.max()))
        imd.set(clim=(dilated.min(), dilated.max()))
        fig.canvas.draw()
        fig.canvas.manager.set_window_title(f"Status: {data['status']}")

    def key_events(event):
        if event.key in keys:
            data["r"] += keys[event.key]
            redraw(ignore_pause=False)
        elif event.key == "enter":
            redraw(ignore_pause=True)
        elif event.key == " ":
            if data["status"] == "pause":
                data["status"] = "waiting"
                fig.canvas.manager.set_window_title(f"Status: {data['status']}")
                redraw()
            else:
                data["status"] = "pause"
                fig.canvas.manager.set_window_title(f"Status: {data['status']}")
        ax.set_title(title % (data["alt"], data["r"]))
        fig.canvas.draw()
    fig.canvas.mpl_connect("key_press_event", key_events)

    def flood_pick(event):
        if event.dblclick:
            data["x"] = np.abs(x-event.xdata+res/2).argmin()
            data["y"] = np.abs(y-event.ydata+res/2).argmin()
            redraw(ignore_pause=True)
    fig.canvas.mpl_connect("button_press_event", flood_pick)
    
    plt.show()
    return (data["x"], data["y"]), float(data["alt"]), data["r"]

def uniform_grid_interp(x, y, Z, xZ=None, yZ=None):
    if xZ is not None:
        x = (x-x[0])/(x[-1]-x[0])*(Z.shape[1]-1)
    if yZ is not None:
        y = (y-y[0])/(y[-1]-y[0])*(Z.shape[0]-1)
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

def gkern(l, sig):
    """creates gaussian kernel with side length `l` and a sigma of `sig`"""
    ax = np.linspace(-(l - 1) / 2., (l - 1) / 2., l)
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

def isotropic_dilation(M, r):
    for _ in range(r):
        M = isodil(M)
    return M

