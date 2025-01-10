from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


outdir = Path(__file__).absolute().parent / "_output"
print(outdir)
outdir = Path("_out")
print(outdir.absolute())

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
x = np.linspace(extent[0], extent[1], num=gridsdata["nx,ny"][0], endpoint=True)
y = np.linspace(extent[0], extent[1], num=gridsdata["nx,ny"][0], endpoint=True)

h, hu, hv, z = np.fromfile(outdir/f"fgout0001.b{1:0>4}", dtype=np.float64).reshape((4, *gridsdata["nx,ny"]), order="F")
fig, axes = plt.subplots(ncols=1, nrows=2, sharex=True, sharey=False, gridspec_kw=dict(hspace=0), layout="tight")
# imh = axes[0, 0].imshow([[0]], extent=extent, origin="lower")
# axes[0,0].set_title("h")
# imhu = axes[0, 1].imshow([[0]], extent=extent, origin="lower")
# axes[0,1].set_title("hu")
# imhv = axes[1, 0].imshow([[0]], extent=extent, origin="lower")
# axes[1,0].set_title("hv")
# imz = axes[1, 1].imshow([[0]], extent=extent, origin="lower")
# axes[1,1].set_title("z")
# imh.set_clim(0, 1)
# imhu.set_clim(-10, 10)
# imhv.set_clim(-10, 10)
# imz.set_clim(0, 1)

# lineh,  = axes[2, 0].plot(x, z.mean(axis=1))
# axes[2,0].set_title("h")
# linehu, = axes[2, 1].plot(x, hu.mean(axis=1))
# linehv, = axes[2, 1].plot(x, hv.mean(axis=1))
# axes[2,1].set_title("flux")
# axes[2, 0].fill_between(x, (z-h).mean(axis=1), color="sienna")

lineh,  = axes[0].plot(x, z.mean(axis=1))
axes[0].set_ylabel("h")
linehu, = axes[1].plot(x, hu.mean(axis=1), label="$hu$")
linehv, = axes[1].plot(x, hv.mean(axis=1), label="$hv$")
axes[1].set_ylabel("flux")
axes[1].set_xlabel("$x$")
axes[0].fill_between(x, (z-h).mean(axis=1), color="sienna")
axes[1].set_ylim(-0.5, 0.5)
axes[1].legend(loc="upper right")

plt.show(block=False)

def read(i):
    fig.suptitle(f"t = {times[i-1]:.2f}")
    h, hu, hv, z = np.fromfile(outdir/f"fgout0001.b{i:0>4}", dtype=np.float64).reshape((4, *gridsdata["nx,ny"]), order="F")
    # imh.set_data(np.ma.MaskedArray(h.T, mask=h.T<1e-5))
    # imhu.set_data(np.ma.MaskedArray(hv.T, mask=h.T<1e-5))
    # imhv.set_data(np.ma.MaskedArray(hu.T, mask=h.T<1e-5))
    # imz.set_data(z.T)
    mid = h.shape[1]//2
    lineh.set_data(x, z.T[mid])
    linehu.set_data(x, hu.T[mid])
    linehv.set_data(x, hv.T[mid])
    # for ax in axes[-1, :]:
    #     ax.relim()
    #     ax.autoscale_view()

anim = FuncAnimation(fig, read, frames=range(1, times.size+1), interval=300)
plt.show()
