from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


outdir = Path(__file__).absolute().parent / "_output5"

def read_clawdata(outdir=outdir):
    clawdata_trans = dict(T=True, F=False)
    clawdata = dict()
    with open(outdir/"claw.data") as file:
        lines = [l for l in file.readlines() if "=" in l]
        for line in lines:
            value, key = line.split("=: ")
            key = key[:key.find(" ")]
            value = [v for v in value.split() if v]
            for e, element in enumerate(value):
                try:
                    value[e] = eval(element)
                except Exception:
                    value[e] = clawdata_trans[element]
            if len(value) == 1:
                value = value[0]
            clawdata[key] = value
    return clawdata

def read_fgrids(outdir=outdir):
    with open(outdir/"fgout_grids.data") as file:
        lines = [l.strip().split("#") for l in file.readlines()[5:] if "#" in l]
        gridsdata = {
            k.strip(): [eval(v) for v in values.split() if v]
            for values, k in lines
        }
        gridsdata = {k: (v[0] if len(v)==1 else v) for k, v in gridsdata.items()}
    return gridsdata

clawdata = read_clawdata()
gridsdata = read_fgrids()
times = np.linspace(0, clawdata["tfinal"], num=clawdata["num_output_times"], endpoint=True)
extent = list(zip(gridsdata["x1, y1"], gridsdata["x2, y2"]))
extent = extent[0] + extent[1]

fig, axes = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True, gridspec_kw=dict(hspace=0))
imh = axes[0, 0].imshow([[0]], extent=extent, origin="lower")
imhu = axes[0, 1].imshow([[0]], extent=extent, origin="lower")
imhv = axes[1, 0].imshow([[0]], extent=extent, origin="lower")
imz = axes[1, 1].imshow([[0]], extent=extent, origin="lower")
imh.set_clim(0, 2)
imhu.set_clim(-10, 10)
imhv.set_clim(-10, 10)
imz.set_clim(1000, 3000)
plt.show(block=False)

def read(i):
    fig.suptitle(f"t = {times[i-1]}")
    h, hu, hv, z = np.fromfile(outdir/f"fgout0001.b{i:0>4}", dtype=np.float64).reshape((4, *gridsdata["nx,ny"]), order="F")
    imh.set_data(np.ma.MaskedArray(h.T, mask=h.T<1e-5))
    imhu.set_data(np.ma.MaskedArray(hv.T, mask=h.T<1e-5))
    imhv.set_data(np.ma.MaskedArray(hu.T, mask=h.T<1e-5))
    imz.set_data(z.T)

anim = FuncAnimation(fig, read, frames=range(1, times.size+1), interval=300)
plt.show()