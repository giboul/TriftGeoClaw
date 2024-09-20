from pathlib import Path
from argparse import ArgumentParser
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from clawpack.pyclaw.solution import Solution


def plot_output(directory, fmt, animate):

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

    # image details
    xmin = solutions[0].states[0].patch.grid.c_centers[0][0, 0]
    ymin = solutions[0].states[0].patch.grid.c_centers[1][0, 0]
    xmax = solutions[0].states[0].patch.grid.c_centers[0][-1, -1]
    ymax = solutions[0].states[0].patch.grid.c_centers[1][-1, -1]
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

    # mutable incremented variable on arrow events
    solix = dict(i=0)
    def increment_solix(event):
        if event.key == "right":
            solix['i'] += 1
        elif event.key == "left":
            solix['i'] -= 1
        solix['i'] = solix['i'] % nsolutions
    fig.canvas.mpl_connect('key_press_event', increment_solix)

    def update(i):
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
                      interpolation="nearest", origin="lower", cmap="Blues", vmin=0, vmax=hmax)
        # Plot empty array in order to show the full extent
        # https://github.com/matplotlib/matplotlib/issues/19374
        ax.imshow(mask, extent=im_extent)
        fig.canvas.draw()
        fig.canvas.flush_events()


    def iplot():
        while plt.fignum_exists(fig.number):
            # Loop over images
            update(solix['i'])
            fig.waitforbuttonpress()

    if animate is False:
        iplot()
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

