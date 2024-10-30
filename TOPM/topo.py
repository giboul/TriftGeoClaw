#!/usr/bin/env python
from pathlib import Path
from argparse import ArgumentParser
from shutil import rmtree
from json import load
from yaml import safe_load
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from skimage.morphology import flood, isotropic_dilation
from skimage.measure import find_contours
from tifffile import TiffFile

projdir = Path().absolute().parent
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    topoconfig = config["TOPM"]

def write_topo(plot=False):

    # Temporary directory
    path = projdir / topoconfig["topography"]
    tempdir = Path("_temp")
    tempdir.mkdir(exist_ok=True)

    print(f"\tINFO: Opening {path}... ")
    with TiffFile(path) as tif:
        Ztif = tif.asarray()
        nx = tif.pages[0].tags["ImageWidth"].value
        ny = tif.pages[0].tags["ImageLength"].value
        x_res, y_res, z_res = tif.pages[0].tags["ModelPixelScaleTag"].value
        xmin = tif.pages[0].tags["ModelTiepointTag"].value[3]
        ymax = tif.pages[0].tags["ModelTiepointTag"].value[4]
    xtif = xmin + (np.arange(nx) + 0.5)*x_res
    ytif = ymax - (np.arange(ny) + 0.5)*y_res

    xmin = topoconfig['bounds']['xmin']
    xmax = topoconfig['bounds']['xmax']
    ymin = topoconfig['bounds']['ymin']
    ymax = topoconfig['bounds']['ymax']
    print(f"\tINFO: Cropping to {xmin, ymin, xmax, ymax = }")
    xmask = (xmin <= xtif) & (xtif <= xmax)
    ymask = (ymin <= ytif) & (ytif <= ymax)
    xtif = xtif[xmask]
    ytif = ytif[ymask][::-1]
    Ztif = Ztif[ymask, :][:, xmask]

    print(f"\tINFO: Downscaling to resolution = {topoconfig['resolution']}")
    x = np.arange(xmin, xmax+topoconfig['resolution'], step=topoconfig['resolution'])
    y = np.arange(ymin, ymax+topoconfig['resolution'], step=topoconfig['resolution'])
    Z = grid_interp(xtif, ytif, Ztif, x, y)
    X, Y = np.meshgrid(x, y)
    mask = dam_mask(X, Y, Z)
    Z[mask] = topoconfig['dam_alt']
    path = projdir / topoconfig["bathymetry"]
    print(f"\tINFO: Saving {path}... ")
    asc_header = "\n".join((
        f"{x.size} ncols",
        f"{y.size} nrows",
        f"{x.min()} xllcenter",
        f"{y.min()} yllcenter",
        f"{topoconfig['resolution']} cellsize",
        f"{999999} nodata_value"
    ))
    np.savetxt(path, Z.flatten(), header=asc_header, comments="")
    print(f"\tINFO: File size is {path.stat().st_size:.2g} bytes.")
    
    print("\tINFO: writing dam coordinates")
    xd = np.linspace(X[mask].min(), X[mask].max(), 100)
    yd = (dam_downstream(xd)+dam_upstream(xd))/2
    mask = np.isfinite(yd)
    np.savetxt("dam.xy", np.vstack((xd[mask], yd[mask])).T)

    rmtree(tempdir)

    y = y[::-1]
    if 'flood_seed' in topoconfig:
        seed = (np.abs(x-topoconfig['flood_seed'][0]).argmin(),
                np.abs(y-topoconfig['flood_seed'][1]).argmin())
    else:
        seed, r = pick_seed(Z, x, y, topoconfig['resolution'])
    flooded = fill_lake(Z, seed[::-1], topoconfig['lake_alt'])
    xc, yc = find_contours(flooded.T, 0.5)[0].T
    xc = xmin + xc/x.size * (xmax - xmin)
    yc = ymin + yc/y.size * (ymax - ymin)
    yc = ymin + (ymax-yc)  # Reverse y
    contour_coords = np.vstack((xc, yc)).T
    np.savetxt(projdir / "TOPM" / "contour.xy", contour_coords)

    dilated = isotropic_dilation(flooded, 50/topoconfig['resolution'])
    xc, yc = find_contours(dilated.T, 0.5)[0].T
    xc = xmin + xc/x.size * (xmax - xmin)
    yc = ymin + yc/y.size * (ymax - ymin)
    yc = ymin + (ymax-yc)  # Reverse y
    dilated_contour_coords = np.vstack((xc, yc)).T
    np.savetxt(projdir / "TOPM" / "contour_dilated.xy", dilated_contour_coords)

    avacs = read_geojson(projdir / topoconfig["avalanches"])
    np.savetxt(projdir / "TOPM" / "avalanches.csv", avacs)

    if plot:
        extent = (xmin, xmax, ymin, ymax) 
        plt.title("Processed topography")
        plt.imshow(Z, extent=extent)
        for i in np.unique(avacs[0]):
            av = avacs.T[i==avacs[0]].T
            plt.fill(*av[1:])
        plt.plot(*contour_coords.T, '-', label="Lake contour")
        plt.plot(*dilated_contour_coords.T, '-', label="Dilated contour")
        plt.scatter(x[seed[0]], y[seed[1]], c="g", label="Fill seed")
        plt.plot(xd, yd, label="Dam middle line")
        plt.plot(xd, dam_upstream(xd), label="Dam upper line")
        plt.plot(xd, dam_downstream(xd), label="Dam lower line")
        plt.legend()
        plt.show()

def expand_bounds(xmin, xmax, ymin, ymax, margin=5*topoconfig['resolution']):
    _xmin = xmin - margin
    _ymin = ymin - margin
    _xmax = xmax + margin
    _ymax = ymax + margin
    return _xmin, _xmax, _ymin, _ymax

def grid_interp(xt, yt, Zt, x, y):
    x = np.clip(x, xt.min(), xt.max())
    y = np.clip(y, yt.min(), yt.max())
    ix = np.clip(np.searchsorted(xt, x)-1, 0, xt.size-2)
    iy = np.clip(np.searchsorted(yt, y)-1, 0, yt.size-2)
    iyz = yt.size-2-iy
    x1 = xt[ix]
    x2 = xt[ix+1]
    y1 = yt[iy]
    y2 = yt[iy+1]
    Z = ((
        + (x2-x)*((y2-y)*Zt[iyz+1, :][:, ix].T).T
        + (x2-x)*((y-y1)*Zt[iyz, :][:, ix].T).T
        + (x-x1)*((y2-y)*Zt[iyz+1, :][:, ix+1].T).T
        + (x-x1)*((y-y1)*Zt[iyz, :][:, ix+1].T).T
    ).T/(y2-y1)).T/(x2-x1)
    return Z[::-1, :]

def dam_mask(x, y, z):
    dam_y1 = dam_upstream(x)
    dam_y2 = dam_downstream(x)
    y = y[::-1]
    return (dam_y1 <= y) & (y <= dam_y2) & (z < topoconfig['dam_alt'])

def dam_upstream(x, offset=0, y0=1171960, x0=2669850, x1=2670561, ymax=1172000):
    x = x + offset
    yd = y0 - 0.3*(x-x0) - 50000/(x-x1) - 300 + offset
    yd[(x < x0) | (x > x1) | (yd > ymax)] = float("inf")
    return yd

def dam_downstream(x, thk=50):
    d = dam_upstream(x)
    u = dam_upstream(x, offset=thk)
    return np.maximum(u, d)

def pick_seed(z_im, x, y, res=0):

    extent = x.min(), x.max(), y.min(), y.max()
    fig, ax = plt.subplots()
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

    data = dict(alt="", x=0, y=0, r=0, status="waiting")
    keys = dict(up=1, right=1, down=-1, left=-1)
    ax.set_title(title % (data["alt"], data["r"]))

    def redraw(ignore_pause=False):
        if data["status"] == "pause" and ignore_pause is False:
            return None
        fig.canvas.manager.set_window_title(f"Flooding...")
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
        if event.key.isnumeric() or event.key == ".":
            data["alt"] += event.key
        elif event.key == "backspace":
            data["alt"] = data["alt"][:-1]
        if event.key in keys:
            data["r"] += keys[event.key]
            redraw(ignore_pause=False)
        elif event.key == "enter":
            redraw(ignore_pause=True)
        elif event.key == " ":
            if data["status"] == "pause":
                data["status"] = "waiting"
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
    return (data["x"], data["y"]), data["r"]

def fill_lake(topo, seed, max_level=0):
    mask = topo < max_level
    initial_value = mask[*seed]
    mask[*seed] = True
    flooded = flood(mask, seed)
    flooded[*seed] = initial_value
    topo[flooded] = max_level
    return flooded

def read_geojson(path):
    with open(path, "r") as file:
        data = load(file)

    coords = []
    for e in data["features"]:
        ix = e['properties']['id']
        points = e['geometry']['coordinates'][0][0]
        for x, y in points:
            coords.append([ix-1, x, y])

    avacs = np.array(coords).T
    return avacs

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-p", "--plot", action="store_true")
    args = parser.parse_args()
    write_topo(plot=args.plot)

