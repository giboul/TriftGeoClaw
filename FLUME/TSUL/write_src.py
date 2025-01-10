from pathlib import Path
from yaml import safe_load
import numpy as np
from matplotlib import pyplot as plt
from clawpack.pyclaw.solution import Solution
from clawpack.visclaw.gridtools import grid_output_2d


projdir = Path(__file__).absolute().parents[1]
with open(projdir / "config.yaml") as file:
    TOPM, AVAC, TSUL = safe_load(file).values()

avid = 5

def read_clawdata(outdir=projdir/"AVAC"/f"_output{avid}"):
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
            clawdata[key] = value
    return clawdata

clawdata = read_clawdata()
x1, y1 = clawdata["lower"]
x2, y2 = clawdata["upper"]
print(x1, x2, y1, y2)
X, Y = np.meshgrid(np.arange(x1, x2, 1.), np.arange(y1, y2, 1.0))
print(X.shape, Y.shape)
frame = Solution(0, path=projdir/"AVAC"/f"_output{avid}", file_format=AVAC["out_format"])
print(frame)
q = grid_output_2d(frame, lambda x: x, X, Y, levels=[1], return_ma = True)
print(q.shape)

