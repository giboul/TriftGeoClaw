from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from clawpack.pyclaw import solution


dir = '_output'
files = list(Path(dir).glob('fort.q*'))
print(f"There are {len(files)} files in {dir}")

# Get first frame
frame = solution.Solution(0, path=dir, file_format='ascii')
# Get bed arrays
x = frame.state.grid.x.centers
s = frame.state.aux[0]

# Initialize figure
fig, ax = plt.subplots(figsize=(5, 2.7), layout="tight")
# Bed and normal depth
water, = ax.plot(x, s, label='Normal depth')
ax.fill_between(x, s, s.min()-.5, color="k", label="Bed")
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$h(x, t)$')

# Function for water line udpating
def update(i):
    frame = solution.Solution(i, path=dir, file_format='ascii')
    h = frame.state.q[0]
    water.set_data(x, h+s)
    water.set_label(rf"$ t = {frame.t:.1f}$ s")
    ax.legend(loc='upper right')
    fig.canvas.draw()
    fig.canvas.flush_events()
# Animate
anim = FuncAnimation(fig, update, frames=len(files), interval=500/len(files))
if False:
    plt.show()
elif True:
    anim.save("movie.gif", fps=len(files)/5, dpi=100)
else:
    del anim
    nstr = len(str(len(files)))
    for i, f in enumerate(files):
        update(i)
        print(f"\rSaving _plots/fig{i:>0{nstr}}.svg...", end="")
        fig.savefig(f"_plots/fig{i:>0{nstr}}.svg")
    print("\nFigures saved")
