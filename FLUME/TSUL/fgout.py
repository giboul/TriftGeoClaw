from pathlib import Path
from yaml import safe_load
from clawpack.geoclaw import fgout_tools
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.path import Path as mPath

projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    config = safe_load(file)["AVAC"]


xc, yc = np.loadtxt(projdir / "TOPM" / "contour1.xy").T
xd, yd = np.loadtxt(projdir / "TOPM" / "contour2.xy").T

fgno = 1
outdir = Path().absolute().parent / "AVAC" / '_output5'
output_format = 'binary'
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, output_format)
fgout_grid.read_fgout_grids_data()

fgframes = list(sorted(outdir.glob("fgout0001.q*")))
outframes = list(sorted(outdir.glob("fort.q*")))
for fgframe in fgframes:
    # asc_header = 
    fgframe = int(str(fgframe)[-4:])
    fgout = fgout_grid.read_frame(fgframe)
    np.savetxt(f"_inflows5/h{fgframe:0>4}.txt", fgout.h.T[::-1])
    plt.figure(1, figsize=(8,8))
    plt.imshow(fgout.B.T[::-1], extent=fgout.extent_edges, cmap=plt.cm.Blues)
    plt.imshow(np.ma.MaskedArray(fgout.eta, fgout.h>0).T[::-1], extent=fgout.extent_edges, cmap=plt.cm.viridis)
    plt.title('Surface at time %s' % fgout.t)
    # fname = f'fgout_frame{fgframe:0>4}.png'
    # plt.savefig(fname)
    plt.show()
