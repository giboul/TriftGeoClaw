from pathlib import Path
from argparse import ArgumentParser
from matplotlib import pyplot as plt
import numpy as np
from params import out_format
from clawpack.pyclaw.solution import Solution


def plot_output(outdir=None):

    # List output files
    outdir = Path(outdir)
    frame_files = list(outdir.glob("fort.q*"))
    nsolutions = len(frame_files)
    solutions = [
        Solution(f, path=outdir, file_format=out_format)
        for f in range(nsolutions)
    ]

    # Initiate figure
    fig, ax = plt.subplots()
    ax.set_xlabel(r"$x$ [m]")
    ax.set_ylabel(r"$y$ [m]")
    title = "t = %.2f"
    ax.set_title(title % 0.)
    ax.set_aspect("equal")
    fig.show()

    # image details
    X, Y = solutions[0].states[0].patch.grid.c_centers
    h0, *_, eta0  = solutions[0].states[0].q
    xmin = X[0, 0]
    ymin = Y[0, 0]
    xmax = X[-1, -1]
    ymax = Y[-1, -1]
    zmin = (eta0-h0).min()
    zmax = (eta0-h0).max()
    hmax = max([state.q[0].max() for sol in solutions for state in sol.states])
    im_extent = xmin, xmax, ymin, ymax
    mask = np.ma.MaskedArray(np.zeros(solutions[0].states[0].q[0].shape), mask=True)

    # mutable incremented variable on arrow events
    solix = dict(i=0)
    def increment_solix(event):
        if event.key == "right":
            solix['i'] += 1
        elif event.key == "left":
            solix['i'] -= 1
        solix['i'] = solix['i'] % nsolutions
    fig.canvas.mpl_connect('key_press_event', increment_solix)

    while plt.fignum_exists(fig.number):
        # Loop over images
        sol = solutions[solix['i']]
        ax.set_title(title % sol.t)
        for artist in ax.collections + ax.lines:
            artist.remove()
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
        plt.waitforbuttonpress()

def read_args():
    parser = ArgumentParser()
    parser.add_argument("outdir", default="_output", nargs="?")
    return parser.parse_args()

if __name__ == "__main__":
    plot_output(read_args().outdir)

