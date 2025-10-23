from argparse import ArgumentParser
from pathlib import Path
import numpy as np
from clawpack.visclaw import gridtools
from clawpack.clawutil.data import ClawData
from clawpack.pyclaw.solution import Solution

from utils import config

def extract(i, x, y, outdir):
    frame_sol = Solution(i, path=outdir, file_format=config["output_format"])
    q = gridtools.grid_output_2d(
        frame_sol,
        lambda q: q,
        x, y,
        levels = "all",
        return_ma=True
    )
    return q, frame_sol.t

def write(AVAC_outdir, extent, outdir, n):

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

    Path(outdir).mkdir(exist_ok=True)
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
    parser = ArgumentParser()
    parser.add_argument("AVAC_outdir", default=config["AVAC_outdir"], type=str, nargs="?")
    parser.add_argument("outdir", default="_output", type=str, nargs="?")
    parser.add_argument("extent", default=None, nargs="?")
    parser.add_argument("n", default=config["bc_size"], type=int, nargs="?")
    return parser.parse_args()

def main():
    args = parse_args()
    write(**args.__dict__)

if __name__ == "__main__":
    main()
