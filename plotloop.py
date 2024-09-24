from pathlib import Path
from argparse import ArgumentParser
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from clawpack.visclaw.Iplotclaw import Iplotclaw
from clawpack.visclaw.frametools import plotframe


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


def plot_output(setplot='setplot.py',
                outdir=None,
                completekey='tab',
                stdin=None,
                stdout=None,
                simple=False,
                controller=None,
                animate=False):
    """Take advantage of Iplotclaw."""
    outdir = Path(outdir or "_output")
    frame_files = list(outdir.glob("fort.q*"))
    nsolutions = len(frame_files)

    # Set up counter
    counter = Counter(max=nsolutions)
    def wait_update_key(event):
        """If keys is in Couter.keys, increment index and update plot."""
        if event.key in counter.keys:
            counter.increment(event.key)
            update(counter.ix)

    ip = Iplotclaw(setplot='setplot.py',
                   outdir=outdir,
                   completekey=completekey,
                   stdin=stdin,
                   stdout=stdout,
                   simple=simple,
                   controller=controller)

    def update(i):
        """Clear axis, show solution `i` and draw."""
        ip.frameno = i
        plotframe(i, ip.plotdata, simple=ip.simple, refresh=True)

    if animate is False:
        plt.close()  # TODO prevent initial empty figure
        update(0)
        fig = plt.gcf()
        fignum = fig.number
        fig.canvas.mpl_connect('key_press_event', wait_update_key)
        while plt.fignum_exists(fignum):
            fig.waitforbuttonpress()
    else:
        anim = FuncAnimation(fig, update, nsolutions, interval=200)
        anim.save("test.mp4")
        plt.show()

def read_args():
    parser = ArgumentParser()
    parser.add_argument("-d", "--outdir", default="_output", nargs="?")
    parser.add_argument("-f", "--fmt", default="_output", nargs="?")
    parser.add_argument("-a", "--animate", action="store_true")
    return parser.parse_args()

if __name__ == "__main__":
    plot_output()

