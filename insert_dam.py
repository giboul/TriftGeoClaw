from pathlib import Path
import numpy as np
import geopandas as gpd
from yaml import safe_load
import topo_utils
from matplotlib.path import Path as mPath
from matplotlib import pyplot as plt


with open("config.yaml") as file:
    config = safe_load(file)

dams = gpd.read_file(Path(config["dam"]).expanduser()).geometry
dam_alts = config["dam_alt"]
if not np.iterable(dam_alts):
    dam_alts = [dam_alts]

x, y, Z = topo_utils.read_asc(config["base_bathymetry"])
topo_shape = Z.shape
X, Y = np.meshgrid(x, y)
Z = Z.flatten()
X = X.flatten()
Y = Y.flatten()
XY = np.column_stack((X, Y))

for dam, dam_alt in zip(dams, dam_alts):
    dam = np.array(dam.exterior.coords)
    in_dam = mPath(dam).contains_points(XY)
    Z[in_dam] = np.maximum(Z[in_dam], dam_alt)

Z = Z.reshape(topo_shape)

topo_utils.write_asc(Path(config["bathymetry"]).expanduser(), Z, x.min(), y.min(), x[1]-x[0])

plt.imshow(Z, extent = (x.min(), x.max(), y.min(), y.max()))
plt.show()
