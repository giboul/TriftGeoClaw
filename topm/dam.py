import numpy as np
from matplotlib.path import Path as mPath
import utils


def dam_mask(pts, X, Y, Z, dam_alt):
    mask = mPath(pts).contains_points(np.column_stack((X.flatten(), Y.flatten()))).reshape(X.shape)
    mask = mask & (Z <= dam_alt)
    return mask


def main():
    from pathlib import Path
    from yaml import safe_load
    from matplotlib import pyplot as plt

    projdir = Path(__file__).parents[1]
    with open(projdir / "config.yaml") as file:
        TOPM = safe_load(file)["TOPM"]
    ix, xd, yd = utils.read_geojson(projdir / "dam.geojson")
    x, y, Z = utils.read_asc(projdir / TOPM["bathymetry"])
    extent = x[0], x[-1], y[-1], y[0]
    x, y = np.meshgrid(x, y)
    mask = dam_mask(np.column_stack((xd, yd)), x, y, Z, TOPM["dam_alt"])
    Z[mask] = TOPM["dam_alt"]
    plt.imshow(Z, extent=extent)
    plt.imshow(np.where(mask, True, np.nan), extent=extent)
    plt.plot(xd, yd)
    plt.show()


if __name__ == "__main__":
    main()

