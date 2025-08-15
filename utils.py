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


def read_times(outdir, gridno):
    outdir = Path(outdir)
    return [
        read_clawdata(p, sep=" ")["time"]
        for p in sorted(outdir.glob(f"fgout{gridno:0>4}.t*"))
    ]


def read_clawdata(path: Path, sep="=: ", comments="#", skiprows=0): # TODO read grids correctly (e.g. with fgno)
    path = Path(path)
    clawdata_trans = dict(T=True, F=False)
    clawdata = dict()
    lines = [line for line in path.read_text().split("\n")[skiprows:] if sep in line]
    for line in lines:
        value, key = [e for e in line.split(sep) if e]
        key = key.strip()
        if comments in key:
            key = key[:key.find(comments)].strip()
        value = [v for v in value.split() if v]
        for e, element in enumerate(value):
            try:
                value[e] = eval(element)
            except Exception:
                value[e] = clawdata_trans.get(element, element)
        clawdata[key] = value[0] if len(value)==1 else value
    return clawdata


def read_datafiles(outdir: Path, ext=".data", *args, **kwargs):
    return {f.stem: read_clawdata(f, *args, **kwargs) for f in outdir.glob(f"*{ext}")}


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

def set_transparent_cmaps(eps, water_color=water_color):

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
        'red': ((0.0, 1.0,1.0,),
                (0.2, 1.0, 1.0),
                (0.5, water_color[0], water_color[0]),
                (1.0, 0.0,0.0,),),
        'green': ((0.0, 0.0,0.0,), 
                (0.2, 0.0, 0.0),
                  (0.5, water_color[1], water_color[1]),
                  (1.0, 0.0,0.0,),),
        'blue': ((0.0, 0.0,0.0,),
                (0.2, 0.0, 0.0),
                 (0.5, water_color[2], water_color[2]),
                 (1.0, 1.0,1.0,),),
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
                  (0.5-eps, 1.0, 1.0),
                  (0.5, 0.0, 0.0),
                  (0.5+eps, 1.0, 1.0),
                  (1.0, 1.0, 1.0),),
    }))

def read_world_image(path):
    path = Path(path)
    im = plt.imread(path)
    ny, nx, _ = im.shape
    world_text = path.with_suffix(".pgw").read_text().strip().split("\n")
    dx, _, _, dy, xul, yul = [float(f) for f in world_text]
    return im, (xul, xul+nx*dx, yul+ny*dy, yul)

