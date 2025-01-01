#!/usr/bin/env python
from pathlib import Path
from argparse import ArgumentParser
from shutil import rmtree
from json import load
from yaml import safe_load
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.widgets import TextBox
from skimage.morphology import flood, isotropic_dilation, isotropic_erosion, disk
from skimage.measure import find_contours
from tifffile import TiffFile

projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    TOPM = config["TOPM"]

def write_topo(plot=False):

    # Temporary directory
    path = projdir / TOPM["topography"]
    tempdir = projdir / "TOPM" / "_temp"
    tempdir.mkdir(exist_ok=True)

    print(f"\tINFO: Opening {path}... ")
    with TiffFile(path) as tif:
        Z = tif.asarray().astype(np.float16)
        nx = tif.pages[0].tags["ImageWidth"].value
        ny = tif.pages[0].tags["ImageLength"].value
        x_res, y_res, z_res = tif.pages[0].tags["ModelPixelScaleTag"].value
        xmin = tif.pages[0].tags["ModelTiepointTag"].value[3]
        ymax = tif.pages[0].tags["ModelTiepointTag"].value[4]
        tifres = tif.pages[0].tags["ModelPixelScaleTag"].value[0]
        # for k, v in tif.pages[0].tags.items():
        #     print(f"\t\t{k}: {v}")
    x = xmin + (np.arange(nx) + 0.5)*x_res
    y = ymax - (np.arange(ny) + 0.5)*y_res

    xmin = TOPM['bounds']['xmin']
    xmax = TOPM['bounds']['xmax']
    ymin = TOPM['bounds']['ymin']
    ymax = TOPM['bounds']['ymax']
    print(f"\tINFO: Cropping to {xmin, ymin, xmax, ymax = }")
    xmask = (xmin <= x) & (x <= xmax)
    ymask = (ymin <= y) & (y <= ymax)
    x = x[xmask]
    y = y[ymask][::-1]
    Z = Z[ymask, :][:, xmask]
    # plt.imshow(Z)
    # plt.show()

    print(f"\tINFO: Downscaling to resolution = {TOPM['resolution']}")
    if np.isclose(TOPM["resolution"], tifres):
        print("\tINFO: Resolution is already as desired.")
        y = y[::-1]
    else:
        x = np.arange(xmin, xmax+TOPM['resolution'], step=TOPM['resolution'])
        y = np.arange(ymin, ymax+TOPM['resolution'], step=TOPM['resolution'])[::-1]
        Z = grid_interp(x[0], y[-1], Z.shape[1], Z.shape[0], x_res, Z, x, y[::-1])
    # plt.imshow(Z)
    # plt.show()
    natpath = projdir / "TOPM" / 'natural_bathymetry.bin'
    print(f"\tINFO: Saving {natpath} (no dam nor smoothing)... ")
    asc_header = "\n".join((
        f"{x.size} ncols",
        f"{y.size} nrows",
        f"{x.min()} xllcenter",
        f"{y.min()} yllcenter",
        f"{TOPM['resolution']} cellsize",
        f"{999999} nodata_value"
    ))
    with open(natpath.with_suffix(".data"), "w") as file:
        file.write(asc_header)
    # Matrix is written in C order anyway... Hence the '.T'
    # https://stackoverflow.com/questions/37905759/save-numpy-array-as-binary-to-read-from-fortran
    Z.reshape(Z.shape, order="F").astype(np.float16).T.tofile(natpath)
    print(f"\tINFO: {natpath} size is {natpath.stat().st_size:.2g} bytes.")

    print(f"\tINFO: Smoothing topo...")
    smooth_radius = max(1, int(TOPM.get("smooth_radius", 5) / TOPM["resolution"]))
    kernel = disk(smooth_radius)
    # kernel = ~isotropic_erosion(disk(smooth_radius), 1)
    kernel = kernel / kernel.sum()
    Z = conv2d(Z, kernel)
    Z = pad_to_shape(Z, (y.size, x.size))
    # plt.imshow(Z)
    # plt.show()

    print("\tINFO: Inserting dam...")
    X, Y = np.meshgrid(x, y)
    X = X.astype(np.float32)
    Y = Y.astype(np.float32)
    Z = Z.astype(np.float16)
    mask = dam_mask(X, Y, Z)
    Z[mask] = TOPM['dam_alt']
    path = projdir / TOPM["bathymetry"]
    # plt.imshow(Z)
    # plt.show()

    print(f"\tINFO: Saving {path}... ")
    asc_header = "\n".join((
        f"{x.size} ncols",
        f"{y.size} nrows",
        f"{x.min()} xllcenter",
        f"{y.min()} yllcenter",
        f"{TOPM['resolution']} cellsize",
        f"{999999} nodata_value"
    ))
    np.savetxt(path, Z.flatten(), header=asc_header, comments="", fmt="%.8e")
    print(f"\tINFO: File size is {path.stat().st_size:.2g} bytes.")
 
    print("\tINFO: writing dam coordinates")
    xdam = np.linspace(X[mask].min(), X[mask].max(), 100)
    ydam = (dam_downstream(xdam)+dam_upstream(xdam))/2
    mask = np.isfinite(ydam)
    np.savetxt(projdir / "TOPM" / "dam.xy", np.column_stack((xdam[mask], ydam[mask])))

    rmtree(tempdir)

    print("\tINFO: Flooding lake")
    y = y[::-1]
    if 'flood_seed' in TOPM:
        seed = (np.abs(x-TOPM['flood_seed'][0]).argmin(),
                np.abs(y[::-1]-TOPM['flood_seed'][1]).argmin())
        r = 3  # TODO
    else:
        seed, TOPM['lake_alt'], r = pick_seed(Z, x, y[::-1], TOPM['resolution'], TOPM['lake_alt'])
    # Fill topo
    Z_lake = Z.copy()
    flooded = fill_lake(Z_lake, seed[::-1], TOPM["lake_alt"])
    dilated = isotropic_dilation(flooded, r)
    Z_lake[dilated] = TOPM["lake_alt"]
    Z_lake[~dilated] = Z.min()
    # plt.imshow(Z_lake)
    # plt.show()

    print("\tINFO: writing qinit")
    # Write qinit
    np.savetxt(projdir/"TSUL"/"qinit.xyz", np.column_stack((
        X.flatten(), Y.flatten(), Z_lake.flatten()
    )), fmt="%.9e")


    print("\tINFO: writing contours", end=" ", flush=True)
    contour1 = contour(isotropic_erosion(fill_lake(Z, seed[::-1], TOPM['lake_alt']).T, 1))
    contour1 = scale_contour(*contour1.T, x.size, y.size, **TOPM['bounds'])
    np.savetxt(projdir / "TOPM" / "contour1.xy", contour1)

    print("2", end=" ", flush=True)
    contour2 = contour(isotropic_dilation(fill_lake(Z, seed[::-1], TOPM['lake_alt']+TOPM["overhang"]).T, 1))
    contour2 = scale_contour(*contour2.T, x.size, y.size, **TOPM['bounds'])
    np.savetxt(projdir / "TOPM" / "contour2.xy", contour2)

    print("3", end=" ", flush=True)
    contour3 = contour(fill_lake(Z, seed[::-1], TOPM['lake_alt']+40).T)
    contour3 = scale_contour(*contour3.T, x.size, y.size, **TOPM['bounds'])
    np.savetxt(projdir / "TOPM" / "contour3.xy", contour3)

    print(f"\tINFO: Saving extents")
    extent = expand_bounds(
        X[flooded].min(), X[flooded].max(),
        Y[flooded].min(), Y[flooded].max(),
        rel_margin=1/20,
        abs_margin=10
    )
    np.savetxt(projdir/"TOPM"/"lake_extent.txt", extent)

    print(f"\tINFO: Saving avalanches.csv")
    avacs = read_geojson(projdir / TOPM["avalanches"])
    np.savetxt(projdir / "TOPM" / "avalanches.csv", avacs)

    if plot:
        plt.title("Processed topography")
        plt.imshow(Z, extent=TOPM['bounds'].values())
        for i in np.unique(avacs[0]):
            av = avacs.T[i==avacs[0]].T
            plt.fill(*av[1:])
        plt.plot(*contour1.T, '-', label="Lake contour")
        plt.plot(*contour2.T, '-', label="Dilated contour")
        plt.plot(*contour3.T, '-', label="Dilated contour")
        plt.scatter(x[seed[0]], y[seed[1]], c="g", label="Fill seed")
        plt.plot(xdam, ydam, label="Dam middle line")
        plt.plot(xdam, dam_upstream(xdam), label="Dam upper line")
        plt.plot(xdam, dam_downstream(xdam), label="Dam lower line")
        plt.legend()
        plt.show(block=False)
        plt.figure()
        plt.imshow(Z_lake, extent=(x.min(), x.max(), y.min(), y.max()))
        plt.scatter(extent[:2], extent[2:], c='r')
        plt.show()

def contour(mask):
    return find_contours(mask, 0.5)[0]

def scale_contour(x, y, nx, ny, xmin, xmax, ymin, ymax):
    x = xmin + x/nx * (xmax - xmin)
    y = ymin + y/ny * (ymax - ymin)
    y = ymin + (ymax-y)  # Reverse y
    return np.column_stack((x, y))

def expand_bounds(x1, x2, y1, y2, rel_margin=1/50, abs_margin=0):
    dx = (x2 - x1) * rel_margin + abs_margin
    dy = (y2 - y1) * rel_margin + abs_margin
    xmin = x1 - dx
    xmax = x2 + dx
    ymin = y1 - dy
    ymax = y2 + dy
    return xmin, xmax, ymin, ymax


def grid_interp(xll, yll, nx, ny, cs, Zt, x, y):
    xur = xll + nx*cs
    yur = yll + ny*cs
    w = np.clip(0, nx-2, ((x-xll)/(xur-xll)*(nx-1)).astype(np.int64))
    s = np.clip(0, ny-2, ((y-yll)/(yur-yll)*(ny-1)).astype(np.int64))
    xw = xll + (xur-xll)*w/(nx-1)
    ys = yll + (yur-yll)*s/(ny-1)
    xe = xw + (xur-xll)/(nx-1)
    yn = ys + (yur-yll)/(ny-1)
    Z = (
        + (x-xw)*((y-ys)*Zt[:, w+1][s+1, :].T).T
        + (x-xw)*((yn-y)*Zt[:, w+1][s  , :].T).T
        + (xe-x)*((y-ys)*Zt[:, w  ][s+1, :].T).T
        + (xe-x)*((yn-y)*Zt[:, w  ][s  , :].T).T
    ) / cs**2
    return Z

def dam_mask(x, y, z, ymin=1.171e6+650, ymax=1172000):
    mask = (dam_upstream(x) <= y)
    mask = mask & (y <= dam_downstream(x))
    mask = mask & (ymin < y) & (y < ymax)
    return (mask & (z < TOPM['dam_alt']))

def dam_upstream(x, offset=0, y0=1171960, x0=2669850, x1=2670561, eps=1e-5):
    x = x + offset
    # np.divide(500, x-x1, out=np.zeros_like(x), where=~np.isclose(x, x1))
    yd = y0 - 0.3*(x-x0) - np.divide(50000, x-x1, out=np.full_like(x, np.inf, dtype=np.float64), where=(x0<x)&(x+1e-8<x1)) - 300 + offset
    # yd[(x < x0) | (x > x1) | (yd > ymax)] = float("inf")
    return yd

def dam_downstream(x, thk=30):
    d = dam_upstream(x)
    u = dam_upstream(x, offset=thk)
    return np.maximum(u, d)

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

def fill_lake(topo, seed, max_level=0):
    mask = topo < max_level
    initial_value = mask[*seed]
    mask[*seed] = True
    flooded = flood(mask, seed)
    flooded[*seed] = initial_value
    # topo[flooded] = max_level
    return flooded

def norm2_kernel(radius: int):
    x, y = np.meshgrid(np.arange(2*radius), np.arange(2*radius))
    k = np.sqrt((x-radius)**2 + (y-radius)**2)
    return k

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

