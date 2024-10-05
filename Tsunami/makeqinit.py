from params import lake_alt, flood_seed
import numpy as np
from matplotlib import pyplot as plt
from skimage.morphology import flood, isotropic_dilation


def write_qinit(filename = "bathy_with_dam.asc"):

    def before_space(str: str):
        return str[:str.index(" ")]
    # Load bathymetry
    with open(filename, "r") as file:
        nx = int(before_space(file.readline()))
        ny = int(before_space(file.readline()))
        xmin = float(before_space(file.readline()))
        ymin = float(before_space(file.readline()))
        res = float(before_space(file.readline()))
        ndv = float(before_space(file.readline()))

    x = xmin + res*np.arange(nx)
    y = ymin + res*np.arange(ny)[::-1]
    z = np.loadtxt(filename, skiprows=6)

    # Fill topo
    z_lake = z.reshape(ny, nx)
    seed_ix, seed_iy, radius = pick_seed(z_lake.copy(), x, y, res)
    flooded = fill_lake(z_lake, (seed_iy, seed_ix), lake_alt)
    dilated = isotropic_dilation(flooded, radius)
    z_lake[dilated] = lake_alt

    # To .xyz format
    z_lake = z_lake.flatten()
    x, y = np.meshgrid(x, y)
    x = x.flatten()
    y = y.flatten()

    # Write qinit
    np.savetxt("qinit.xyz", np.vstack((x, y, z_lake)).T)


def pick_seed(z_im, x, y, res=0):

    extent = x.min(), x.max(), y.min(), y.max()
    fig, ax = plt.subplots()
    ax.imshow(z_im, extent=extent, cmap="inferno")
    imd = ax.imshow(np.ma.MaskedArray([[1]], mask=True), extent=extent, cmap="Reds")
    imf = ax.imshow(np.ma.MaskedArray([[1]], mask=True), extent=extent, cmap="Blues_r")
    title = "Max altitude: %s   Dilation radius: %i"

    data = dict(alt="", x=0, y=0, r=0)
    keys = dict(up=1, right=1, down=-1, left=-1)
    ax.set_title(title % (data["alt"], data["r"]))

    def redraw():
        ax.set_title(title % (data["alt"], data["r"]))
        flooded = fill_lake(z_im.copy(), (data["y"], data["x"]), float(data["alt"] or 0))
        dilated = isotropic_dilation(flooded, data["r"])
        flooded = np.ma.MaskedArray(z_im, ~flooded)
        dilated = np.ma.MaskedArray(z_im, ~dilated)
        imf.set_data(flooded)
        imd.set_data(dilated)
        imf.set(clim=(flooded.min(), flooded.max()))
        imd.set(clim=(dilated.min(), dilated.max()))
        fig.canvas.draw()

    def enter_alt(event):
        if event.key.isnumeric() or event.key == ".":
            data["alt"] += event.key
            redraw()
        elif event.key == "backspace":
            data["alt"] = data["alt"][:-1]
            redraw()
    fig.canvas.mpl_connect("key_press_event", enter_alt)

    def flood_pick(event):
        if event.dblclick:
            data["x"] = np.abs(x-event.xdata+res/2).argmin()
            data["y"] = np.abs(y-event.ydata+res/2).argmin()
            redraw()
    fig.canvas.mpl_connect("button_press_event", flood_pick)

    def dilate(event):
        if event.key in keys:
            data["r"] += keys[event.key]
            redraw()
    fig.canvas.mpl_connect("key_press_event", dilate)

    plt.show()
    return data["x"], data["y"], data["r"]


def fill_lake(topo, seed, max_level=0):
    mask = topo < max_level
    mask[*seed] = True
    flooded = flood(mask, seed)
    topo[flooded] = max_level
    return flooded


if __name__ == "__main__":
    write_qinit()
