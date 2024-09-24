from pathlib import Path
from argparse import ArgumentParser
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from clawpack.pyclaw.solution import Solution


class Counter:
    """Mutable object to keep track of current plot index."""

    def __init__(self, ix: int = 0, max: int = np.inf):
        self.ix = ix
        self.max = max
        self.keys = dict(right=1, left=-1, up=1, down=-1)

    def increment(self, key):
        """Get increment value (default is 0) and update counter."""
        val = self.keys[key]
        self.ix = (self.ix + val) % self.max


def plot_output(directory, fmt, animate):
    """Create new figure, plot initial solution and set up
    interactive mode or animation."""
    # List output files
    outdir = Path(directory)
    frame_files = list(outdir.glob("fort.q*"))
    nsolutions = len(frame_files)
    solutions = [
        Solution(f, path=outdir, file_format=fmt)
        for f in range(nsolutions)
    ]

    # Initiate figure
    fig, ax = plt.subplots()
    title = "t = %.2f"
    ax.set_title(title % 0.)
    ax.set_aspect("equal")
    fig.show()

    # Get image extents
    xmin = solutions[0].states[0].patch.grid.c_nodes[0][0, 0]
    ymin = solutions[0].states[0].patch.grid.c_nodes[1][0, 0]
    xmax = solutions[0].states[0].patch.grid.c_nodes[0][-1, -1]
    ymax = solutions[0].states[0].patch.grid.c_nodes[1][-1, -1]
    q0 = solutions[0].states[0].q
    for state in solutions[0].states:
        X, Y = state.patch.grid.c_centers
        xmin = min(xmin, X[0, 0])
        ymin = min(ymin, Y[0, 0])
        xmax = max(xmax, X[-1, -1])
        ymax = max(ymax, Y[-1, -1])
    h0, *_, eta0  = solutions[0].states[0].q
    zmin = min((st.q[-1]-st.q[0]).min() for st in solutions[0].states)
    zmax = max((st.q[-1]-st.q[0]).max() for st in solutions[0].states)
    hmax = max(st.q[0].max() for sol in solutions for st in sol.states)
    im_extent = xmin, xmax, ymin, ymax
    mask = np.ma.MaskedArray(np.zeros((1,1)), mask=True)

    # Set up counter
    counter = Counter(max=nsolutions)
    def wait_update_key(event):
        """If keys is in Couter.keys, increment index and update plot."""
        if event.key in counter.keys:
            counter.increment(event.key)
            update(counter.ix)
    fig.canvas.mpl_connect('key_press_event', wait_update_key)

    def update(i):
        """Clear axis, show solution `i` and draw."""
        sol = solutions[i]
        ax.clear()
        ax.set_xlabel(r"$x$ [m]")
        ax.set_ylabel(r"$y$ [m]")
        ax.set_title(title % sol.t)

        for state in sol.states:
            # read solution
            h, hu, hv, *_, eta = state.q
            # adapt extents
            x, y = state.patch.grid.c_nodes
            patch_extent = x[0, 0], x[-1, -1], y[0, 0], y[-1, -1]
            # show
            ax.imshow((eta-h).T, extent=patch_extent, interpolation="nearest",
                      origin="lower", cmap="inferno", vmin=zmin, vmax=zmax)
            ax.imshow(np.ma.MaskedArray(h, mask=h<=1e-3).T, extent=patch_extent,
                      interpolation="nearest", origin="lower",
                      cmap="Blues", vmin=0, vmax=hmax)
        # Plot empty array in order to show the full extent
        # https://github.com/matplotlib/matplotlib/issues/19374
        ax.imshow(mask, extent=im_extent)
        fig.canvas.draw()
        fig.canvas.flush_events()

    if animate is False:
        update(0)
        while plt.fignum_exists(fig.number):
            fig.waitforbuttonpress()
    else:
        anim = FuncAnimation(fig, update, nsolutions, interval=200)
        anim.save("test.mp4")
        plt.show()

def read_args():
    parser = ArgumentParser()
    parser.add_argument("-d", "--directory", default="_output", nargs="?")
    parser.add_argument("-f", "--fmt", default="_output", nargs="?")
    parser.add_argument("-a", "--animate", action="store_true")
    return parser.parse_args()

if __name__ == "__main__":
    plot_output(**read_args().__dict__)

