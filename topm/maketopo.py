from pathlib import Path
import utils
from yaml import safe_load
import numpy as np
from matplotlib import pyplot as plt
from skimage.morphology import flood


projdir = Path(__file__).parent.parent
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    TOPM = config["TOPM"]

print(f"Reading {projdir / TOPM['topography']}")
xt, yt, Z = utils.read_tiff(projdir / TOPM["topography"])

print(f"Print cropping to {TOPM['bounds']}")
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

print(f"Rescaling to {TOPM['resolution']}")
x = np.arange(xmin, xmax+TOPM["resolution"], step=TOPM["resolution"])
y = np.arange(ymin, ymax+TOPM["resolution"], step=TOPM["resolution"])
Z = utils.uniform_grid_interp(x, y, Z, xt, yt)

natpath = projdir / "topm" / 'natural_bathymetry.bin'
print(f"Saving {natpath}")
asc_header = f"""
{x.size} ncols
{y.size} nrows
{x.min()} xllcenter
{y.min()} yllcenter
{TOPM['resolution']} cellsize
{999999} nodata_value
""".strip()
with open(natpath.with_suffix(".data"), "w") as file:
    file.write(asc_header)
Z.reshape(Z.shape, order="F").astype(np.float16).T.tofile(natpath)

# print(f"Smoothing topography")
# smooth_radius = max(1, int(TOPM.get("smooth_radius", 5) / TOPM["resolution"]))
# Z = utils.apply_gkernel(Z, smooth_radius, sig=1)

print("Inserting dam")
def dam(x, y0=1171660, x0=2669850, x1=2670561, eps=1e-5):
    return y0 - 0.3*(x-x0) - np.divide(
            50000, x-x1,
            out=np.full_like(x, np.inf, dtype=np.float64),
            where=((x0 + eps < x) & (x + eps < 2670500))
    )
def dam_mask(x, y, z, y0=1171960, x0=2669850, x1=2670561, eps=1e-5):
    dam_dw = dam(x[0])
    dam_up = dam(x[0]+30)+30
    dam_dm = np.isfinite(dam_dw)
    dam_upm = np.isfinite(dam_up)
    return (dam_dw <= y) & (y <= dam_up) & (z < TOPM['dam_alt'])
y = y[::-1]
X, Y = np.meshgrid(x, y, copy=False)
mask = dam_mask(X, Y, Z)
Z[mask] = TOPM['dam_alt']

print(f"Saving {projdir/'topm'/'bathymetry.asc'}")
bathypath = projdir / "topm" / "bathymetry.asc"
np.savetxt(bathypath, Z.flatten(), header=asc_header, comments="", fmt="%.8e")

print("Filling lake")
if 'flood_seed' in TOPM:
    seed = (np.abs(x-TOPM['flood_seed'][0]).argmin(),
            np.abs(y-TOPM['flood_seed'][1]).argmin())
    r = int(TOPM["dilation_radius"]/TOPM["resolution"])
else:
    seed, TOPM['lake_alt'], r = utils.pick_seed(Z, x, y, TOPM['resolution'], TOPM['lake_alt'])
flooded = utils.fill_lake(Z, seed[::-1], TOPM["lake_alt"])
dilated = utils.isotropic_dilation(flooded, r)
print(f"Saving {projdir / 'tsul' / 'qinit.xyz'}")
Z[flooded] = TOPM["lake_alt"]
np.savetxt(projdir / "tsul" / "qinit.xyz", np.column_stack((
    X.flatten(), Y.flatten(), np.where(dilated, TOPM["lake_alt"], Z.min()).flatten()
)), fmt="%.9e")

contour = utils.find_contour(dilated.T, extent, x.size, y.size)
np.savetxt(projdir / "tsul" / "contour.xy", contour)

extent = utils.expand_bounds(
    X[flooded].min(), X[flooded].max(),
    Y[flooded].min(), Y[flooded].max(),
    rel_margin=1/20,
    abs_margin=10
)
np.savetxt(projdir / "topm" / "lake_extent.txt", extent)

avacs = utils.read_geojson(projdir / TOPM["avalanches"])
np.savetxt(projdir / "avac" / "avalanches.csv", avacs)

if False:
    plt.title("Processed topography")
    plt.imshow(Z, extent=TOPM['bounds'].values())
    for i in np.unique(avacs[0]):
        av = avacs.T[i==avacs[0]].T
        plt.fill(*av[1:])
    plt.plot(*contour.T, '-', label="Lake contour")
    plt.scatter(x[seed[0]], y[seed[1]], c="g", label="Fill seed")
    plt.plot(
        (*extent[:2], *extent[:2][::-1], extent[0]),
        (extent[2], extent[2], extent[3], extent[3], extent[2])
    , c='k', label="TSUL box")
    plt.legend()
    plt.show()
