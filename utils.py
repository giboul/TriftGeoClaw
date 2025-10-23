from pathlib import Path
import json
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from yaml import safe_load


with open("config.yaml") as file:
    config = safe_load(file)


def read_fgout_bin(outdir, nx, ny, frameno, gridno=1, nvars=4):
    outdir = Path(outdir)
    q = np.fromfile(outdir / f"fgout{gridno:0>4}.b{frameno+1:0>4}", np.float64)
    q = q.reshape(nvars, nx, ny, order="F")
    return  q



def read_geojson(path):
    with open(path, "r") as file:
        data = json.load(file)

    coords = []
    for e in data["features"]:
        ix = e['properties']['id']
        points = e['geometry']['coordinates'][0][0]
        for x, y in points:
            coords.append([ix, x, y])

    avacs = np.array(coords).T
    return avacs


water_color = (0.372, 0.588, 0.6)

def set_transparent_cmaps(eps=1e-1, water_color=water_color):

    mpl.colormaps.register(mpl.colors.LinearSegmentedColormap("Reds_water", {
        'red': ((0.0, water_color[0], water_color[0]),
                (0.5, 1.0, 1.0,),
                (1.0, 1.0, 1.0,),),
        'green': ((0.0, water_color[1], water_color[1]),
                  (0.5, 1.0, 1.0,),
                  (1.0, 0.0, 0.0,),),
        'blue': ((0.0, water_color[2], water_color[2]),
                 (0.5, 1.0, 1.0,),
                 (1.0, 0.0, 0.0,),),
    }))


    mpl.colormaps.register(mpl.colors.LinearSegmentedColormap("RdBu_water", {
        'red': ((0.0, 0.0, 0.0,),
                (0.5, water_color[0], water_color[0]),
                (1.0, 1.0, 1.0,),),
        'green': ((0.0, 0.0, 0.0,), 
                  (0.5, water_color[1], water_color[1]),
                  (1.0, 0.0, 0.0,),),
        'blue': ((0.0, 1.0, 1.0,),
                 (0.5, water_color[2], water_color[2]),
                 (1.0, 0.0, 0.0,),),
    }))

    mpl.colormaps.register(mpl.colors.LinearSegmentedColormap("RdBu_tc", {
        'red': ((0.0, 1.0, 1.0),
                (0.5, 1.0, 1.0),
                (1.0, 0.0, 0.0),),
        'green': ((0.0, 0.0, 0.0), 
                  (0.5, 1.0, 1.0),
                  (1.0, 0.0, 0.0),),
        'blue': ((0.0, 0.0, 0.0),
                 (0.5, 1.0, 1.0),
                 (1.0, 1.0, 1.0),),
        'alpha': ((0.0, 1.0, 1.0),
                  (0.5-eps, 1.0, 0.0),
                  (0.5+eps, 0.0, 1.0),
                  (1.0, 1.0, 1.0),),
    }))

def read_world_image(path):
    path = Path(path)
    im = plt.imread(path)
    ny, nx, _ = im.shape
    world_text = path.with_suffix(".pgw").read_text().strip().split("\n")
    dx, _, _, dy, xul, yul = [float(f) for f in world_text]
    return im, (xul, xul+nx*dx, yul+ny*dy, yul)

