from sys import argv
from subprocess import call
from pathlib import Path
from yaml import safe_load

projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    config = safe_load(file)["TSUL"]

if config["inflow"] == "bc":
    module = dict(bc="bc_inflows.py", src="src_inflows.py")[config["inflow"]]
    cmd = ["python", module, *argv[1:]]
    call(cmd)

