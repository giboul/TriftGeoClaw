from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
from pathlib import Path


mx = 1920
my = 1280

def norm(a):
    amin = a.min()
    amax = a.max()
    if np.isclose(amin, amax):
        return np.zeros_like(a)
    return (a-amin)/(amax-amin)
inflow_dir = Path(__file__).parent/"_inflows5"
files = list(inflow_dir.glob("cut*.txt"))
x, y = np.loadtxt(Path(__file__).parents[1]/"TOPM"/"contour2.xy").T

fig, axes = plt.subplots(ncols=3, sharey=True, sharex=True)
axes[0].set_title("$h$")
axes[1].set_title("$hu$")
axes[2].set_title("$hv$")
pts = []
lines = []
for ax in axes:
    ax.set_aspect("equal")
    line = LineCollection(np.c_[x[:-1], y[:-1], x[1:], y[1:]].reshape(-1, 2, 2))
    line.set(linewidths=4)
    line.set(cmap=plt.cm.viridis)
    line.set(array=norm(x))
    ax.add_collection(line)
    ax.update_datalim(np.column_stack((x, y)))
    ax.autoscale_view()
    lines.append(line)

state = dict(i=0)
def update(event):
    global i
    if event.key in ["up","right"]:
        state["i"] += 1
    elif event.key in ["down", "left"]:
        state["i"] -=1
    state["i"] = min(max(0,state["i"]), len(files)-1)
    data = np.loadtxt(inflow_dir/f"cut{state['i']:0>4}.txt")
    for i, line in enumerate(lines):
        line.set_array(norm(data[:,i+2]))
    fig.canvas.draw()
plt.gcf().canvas.mpl_connect("key_press_event", update)
plt.show()
