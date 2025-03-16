from pathlib import Path
from yaml import safe_load
import numpy as np
from matplotlib import pyplot as plt
import utils
from dam import dam_mask


projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    TOPM = config["TOPM"]

bathymetry_path = projdir / TOPM["original_bathymetry"]
if bathymetry_path.suffix in [".tif", ".tiff"]:
    xt, yt, Z = utils.read_tiff(bathymetry_path)
else:
    xt, yt, Z = utils.read_asc(bathymetry_path)

print(f"Print cropping to {TOPM['bounds']}", flush=True)
xmin = TOPM["bounds"]["xmin"]
xmax = TOPM["bounds"]["xmax"]
ymin = TOPM["bounds"]["ymin"]
ymax = TOPM["bounds"]["ymax"]
extent = TOPM["bounds"].values()
xmask = (xmin <= xt) & (xt <= xmax)
ymask = (ymin <= yt) & (yt <= ymax)
xt = xt[xmask]
yt = yt[ymask]
Z = Z[ymask, :][:, xmask]

print(f"Rescaling to {TOPM['resolution']}", flush=True)
x = np.arange(xmin, xmax+TOPM["resolution"], step=TOPM["resolution"])
y = np.arange(ymin, ymax+TOPM["resolution"], step=TOPM["resolution"])
Z = utils.uniform_grid_interp(x, y, Z, xt, yt)

# print(f"Smoothing topography")
# smooth_radius = max(1, int(TOPM.get("smooth_radius", 5) / TOPM["resolution"]))
# Z = utils.apply_gkernel(Z, smooth_radius, sig=1)

dam_path = TOPM.get("dam")
if dam_path is not None:
    dam_path = projdir / dam_path
    print("Inserting dam", flush=True)
    X, Y = np.meshgrid(x, y[::-1], copy=False)
    if dam_path.suffix == ".geojson":
        _, xd, yd = utils.read_geojson(dam_path)
    elif dam_path.suffix == ".csv":
        xd, yd = np.loadtxt(dam_path, delimiter=",").T
    else:
        xd, yd = np.loadtxt(dam_path).T
    mask = dam_mask(np.column_stack((xd, yd)), X, Y, Z, TOPM["dam_alt"])
    Z[mask] = TOPM['dam_alt']

utils.write_asc(projdir/TOPM["bathymetry"],
                Z, x[0], y[0], TOPM["resolution"],
                fmt="%.8e")

if True:
    plt.title("Processed topography")
    plt.imshow(Z, extent=TOPM['bounds'].values())
    plt.show()

