"""
Script for reading results around the boundaries of the wave impulse wave simulation.
The data is written in `outdir`/_bc_inflows with one file per boundary per timestep.
"""
from pathlib import Path
import numpy as np
from numpy.typing import ArrayLike
from clawpack.visclaw import gridtools
from clawpack.clawutil.data import ClawData
from clawpack.pyclaw.solution import Solution

def extract(i: int,
            x: ArrayLike,
            y: ArrayLike,
            outdir: Path="_output",
            output_format: str="ascii"):
    """
    Read the result at timestep `i` accross all boundaries.
    
    Parameters
    ----------
        i: int
            Timestep index (from results)
        x: ArrayLike
            The x-coordinates of the the boundaries
        y: ArrayLike
            The y-coordinates of the boundaries

    Return
    ------
    q: np.ndarray
        the solution matrix
    t: np.ndarray
        the time of the solution
    """
    frame_sol = Solution(i, path=outdir, file_format=output_format)
    q = gridtools.grid_output_2d(
        frame_sol,
        lambda q: q,
        x, y,
        levels = "all",
        return_ma=True
    )
    return q, frame_sol.t

def write(AVAC_outdir, extent, outdir, n):
    """
    Read the results in `AVAC_outdir` and write the 
    results along the boundaries defined by `extent` to `outdir` 
    taking `n` samples per boundary.
    """
    outdir = Path(outdir)
    AVAC_outdir = Path(AVAC_outdir)
    ntimes = len(list(AVAC_outdir.expanduser().glob("fort.q*")))

    boundaries = ("bottom", "right", "top", "left")
    if extent is None:
        clawdata = ClawData()
        clawdata.read("claw.data", force=True)
        xmin, ymin = clawdata["lower"]
        xmax, ymax = clawdata["upper"]
    else:
        xmin, xmax, ymin, ymax = extent
    x = np.hstack((
        np.linspace(xmin, xmax, n, endpoint=True),  # South
        np.full(n, xmax),  # East
        np.linspace(xmax, xmin, n, endpoint=True),  # North
        np.full(n, xmin),  # West
    ))
    y = np.hstack((
        np.full(n, ymin),
        np.linspace(ymin, ymax, n),
        np.full(n, ymax),
        np.linspace(ymax, ymin, n)
    ))

    dist = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
    dist = np.hstack((0, dist, 2*dist[-1]-dist[-2]))

    outdir.mkdir(exist_ok=True)
    for b in boundaries:
        for f in outdir.glob(f"{b}*.npy"):
            f.unlink()

    inflows_dir = outdir / "_bc_inflows"
    inflows_dir.mkdir(exist_ok=True)

    times = []
    for ti in range(ntimes):
        print(f"Saving cut {ti+1:>{4}}/{ntimes}...", end="\r")
        q, t = extract(ti, x, y, AVAC_outdir)
        times.append(t)
        h, hu, hv, eta = np.array(q)
        for bi, boundary in enumerate(boundaries):
            s = slice(bi*n, (bi+1)*n)
            data = np.column_stack((x[s], y[s], h[s], hu[s], hv[s]))
            path = inflows_dir / f"{boundary}_{ti:0>{4}}.npy"
            np.savetxt(path, data, comments="")
            data.reshape(data.shape, order="F").astype(np.float64).T.tofile(path)
    print()

    np.savetxt(inflows_dir / "times.txt", times)


def parse_args():
    from config import config
    from argparse import ArgumentParser
    parser = ArgumentParser()
    config_extent = config["lower"][0], config["upper"][0], config["lower"][1], config["upper"][1]
    parser.add_argument("--extent", "-e", type=float, nargs=4, default=config_extent)
    parser.add_argument("--AVAC_outdir", "-a", type=str, default=config["AVAC_outdir"])
    parser.add_argument("--outdir", "-o", type=str, default="_output")
    parser.add_argument("--output_format", "-f", type=str, default=config["output_format"])
    parser.add_argument("--bc_size", "-n", type=int, default=config["bc_size"])
    return parser.parse_args()


def main():
    args = parse_args()
    write(**args.__dict__)


if __name__ == "__main__":
    main()
