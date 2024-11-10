from sys import argv
from subprocess import call
from pathlib import Path
from yaml import safe_load

projdir = Path().absolute().parent
with open(projdir / "config.yaml") as file:
    config = safe_load(file)["TSUL"]

module = dict(bc="bc_inflows.py", src="src_inflows.py")[config["inflow"]]

cmd = ["python", module, *argv[1:]]
print(f"{cmd = }")
call(cmd)
