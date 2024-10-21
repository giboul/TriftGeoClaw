from params import lake_alt, flood_seed
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
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
    imb = ax.imshow(z_im, extent=extent)
    imd = ax.imshow(np.ma.MaskedArray([[1]], mask=True), extent=extent, cmap="Reds")
    imf = ax.imshow(np.ma.MaskedArray([[1]], mask=True), extent=extent, cmap="Blues")
    title = "Max altitude: %s   Dilation radius: %i"
    print(imb.get_cmap())
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
        flooded = fill_lake(z_im.copy(), (data["y"], data["x"]), float(data["alt"] or 0))
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
    return data["x"], data["y"], data["r"]


def fill_lake(topo, seed, max_level=0):
    mask = topo < max_level
    initial_value = mask[*seed]
    mask[*seed] = True
    flooded = flood(mask, seed)
    flooded[*seed] = initial_value
    topo[flooded] = max_level
    return flooded


if __name__ == "__main__":
    write_qinit()
